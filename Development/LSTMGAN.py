import os
import sys
import numpy as np
import pandas as pd
import tensorflow as tf

TRAIN_DATA_FILE = './PS_Helix_500.csv'					# Dataset location
NUM_EPOCHS = 2000										# Number of training epochs
MAX_ATOMS = 150											# Maximum protein chain length
SEQ_LEN = MAX_ATOMS * 2									# Total prediction sequence length (2 angles per atom)
N_STRUCTURES = 1										# Number of structures to predict

class ModelConfig(object):
	def __init__(self):
		self.num_layers = 1								# Number of LSTM layers
		self.rnn_size = 64								# Number of LSTM units
		self.hidden_size = 32							# Number of hidden layer units
		self.num_mixtures = 4
		self.batch_size = 4
		self.num_steps = 10
		self.dropout_rate = 0.5							# Dropout rate
		self.learning_rate = 0.001						# Learning rate

def load_training_data():
	''' Returns a matrix of training data '''
	data = pd.read_csv(TRAIN_DATA_FILE, index_col=0, sep=';')
	data.drop(data.columns[0], axis=1, inplace=True)	# Remove names
	data = data.div(data.abs().max(axis=1), axis=0)		# Normalize by row
	return data.values

class DataLoader(object):
	def __init__(self, data, batch_size=128, num_steps=1):
		self.batch_size = batch_size
		self.n_data, self.seq_len = data.shape
		self._data = data[:self.batch_size, :]
		self.num_steps = num_steps
		self._data = self._data.reshape((self.batch_size, self.seq_len, 1))
		self._reset_pointer()
	def _reset_pointer(self):
		self.pointer = 0
	def reset(self):
		self._reset_pointer()
	def has_next(self):
		return self.pointer + self.num_steps < self.seq_len - 1
	def next_batch(self):
		batch_xs = self._data[:, self.pointer:self.pointer + self.num_steps, :]
		batch_ys = self._data[:, self.pointer + 1:self.pointer + self.num_steps + 1, :]
		self.pointer = self.pointer + self.num_steps
		return batch_xs, batch_ys

def reset_session_and_model():
	''' Resets the TensorFlow default graph and session '''
	tf.reset_default_graph()
	sess = tf.get_default_session()
	if sess: sess.close()

class MDNModel(object):
	def __init__(self, config, is_training=True):
		self.batch_size = config.batch_size
		self._config = config
		self.rnn_size = config.rnn_size
		self.num_layers = config.num_layers
		self.hidden_size = config.hidden_size
		self.num_steps = config.num_steps
		self.num_mixtures = config.num_mixtures
		self.n_gmm_params = self.num_mixtures * 3
		self.learning_rate = config.learning_rate
		self.is_training = is_training
		with tf.variable_scope('mdn_model', reuse=(not self.is_training)): self._build_model()
	def _build_model(self):
		''' Build the MDN Model '''
		self.x_holder = tf.placeholder(tf.float32, [self.batch_size, self.num_steps, 1], name='x')
		self.y_holder = tf.placeholder(tf.float32, [self.batch_size, self.num_steps, 1], name='y')
		multi_rnn_cell = tf.nn.rnn_cell.MultiRNNCell([tf.nn.rnn_cell.LSTMCell(self.rnn_size) for _ in range(self.num_layers)], state_is_tuple=True)
		self.init_state = multi_rnn_cell.zero_state(self.batch_size, tf.float32)
		rnn_outputs, self.final_state = tf.nn.dynamic_rnn(cell=multi_rnn_cell, inputs=self.x_holder, initial_state=self.init_state)
		w1 = tf.get_variable('w1', shape=[self.rnn_size, self.hidden_size], dtype=tf.float32, initializer=tf.truncated_normal_initializer(stddev=0.2))
		b1 = tf.get_variable('b1', shape=[self.hidden_size], dtype=tf.float32, initializer=tf.constant_initializer())
		h1 = tf.nn.sigmoid(tf.matmul(tf.reshape(rnn_outputs, [-1, self.rnn_size]), w1) + b1)
		w2 = tf.get_variable('w2', shape=[self.hidden_size, self.n_gmm_params], dtype=tf.float32, initializer=tf.truncated_normal_initializer(stddev=0.2))
		b2 = tf.get_variable('b2', shape=[self.n_gmm_params], dtype=tf.float32, initializer=tf.constant_initializer())
		gmm_params = tf.matmul(h1, w2) + b2
		#print(gmm_params)
		mu_ = gmm_params[:, : self.num_mixtures]
		sigma_ = gmm_params[:, self.num_mixtures: 2 * self.num_mixtures]
		pi_ = gmm_params[:, 2 * self.num_mixtures:]
		self.mu = mu_
		self.sigma = tf.exp(sigma_ / 2.0)
		self.pi = tf.nn.softmax(pi_)
		#print(self.mu)
		#print(self.sigma)
		if self.is_training:
			self.optimizer = tf.train.AdamOptimizer(self.learning_rate)
			#print(self.y_holder)
			mixture_p = tf.contrib.distributions.Normal(self.mu, self.sigma).prob(tf.reshape(self.y_holder, (-1, 1)))
			mixture_p = tf.multiply(self.pi, mixture_p)
			output_p = tf.reduce_sum(mixture_p, reduction_indices=1, keepdims=True)
			log_output_p = tf.log(output_p)
			mean_log_output_p = tf.reduce_mean(log_output_p)
			self.loss = -mean_log_output_p
			self.train_op = self.optimizer.minimize(self.loss)
	def train_for_epoch(self, sess, data_loader):
		assert self.is_training, 'Must be training model'
		cur_state = sess.run(self.init_state)
		data_loader.reset()
		epoch_loss = []
		while data_loader.has_next():
			batch_xs, batch_ys = data_loader.next_batch()
			batch_xs = batch_xs.reshape((self.batch_size, self.num_steps, 1))
			batch_ys = batch_ys.reshape((self.batch_size, self.num_steps, 1))
			_, batch_loss_, new_state_ = sess.run([self.train_op, self.loss, self.final_state], feed_dict={self.x_holder: batch_xs, self.y_holder: batch_ys, self.init_state: cur_state,})
			cur_state = new_state_
			epoch_loss.append(batch_loss_)
		return np.mean(epoch_loss)
	def predict(self, sess, seq_len=1000):
		assert not self.is_training, 'Must be testing model'
		cur_state = sess.run(self.init_state)
		preds = []
		preds.append(np.random.uniform())
		for step in range(seq_len):
			batch_xs = np.array(preds[-1]).reshape((self.batch_size, self.num_steps, 1))
			mu_, sigma_, pi_, new_state_ = sess.run([self.mu, self.sigma, self.pi, self.final_state], feed_dict={self.x_holder: batch_xs, self.init_state: cur_state})
			select_mixture = np.random.choice(self.num_mixtures, p=pi_[0])
			new_pred_ = np.random.normal(loc=mu_[0][select_mixture], scale=sigma_[0][select_mixture])
			preds.append(new_pred_)
			cur_state = new_state_
		return preds[1:]

def Run(choice):
	train_config = ModelConfig()
	train_config.learning_rate = 0.0003
	test_config = ModelConfig()
	test_config.batch_size = 1
	test_config.num_steps = 1
	# Train
	if choice == 'train':
		os.makedirs('weights', )
		data = load_training_data()
		reset_session_and_model()
		with tf.Session() as sess:
			loader = DataLoader(data=data, batch_size=train_config.batch_size, num_steps=train_config.num_steps)
			train_model = MDNModel(train_config, True)
			test_model = MDNModel(test_config, False)
			sess.run(tf.global_variables_initializer())
			saver = tf.train.Saver(max_to_keep=0)
			for idx in range(NUM_EPOCHS):
				epoch_loss = train_model.train_for_epoch(sess, loader)
				print('Epoch: {}\tLoss: {}'.format(idx+1, epoch_loss))
			saver.save(sess, f'./weights/mdnmodel.ckpt')
		true_data = data[0]
	# Predict
	elif choice == 'predict':
		ckpt_path = f'./weights/mdnmodel.ckpt'
		for structure in range(N_STRUCTURES):
			reset_session_and_model()
			with tf.Session() as sess:
				test_model = MDNModel(test_config, True)
				test_model.is_training = False
				sess.run(tf.global_variables_initializer())
				saver = tf.train.Saver()
				saver.restore(sess, ckpt_path)
				fake_data = test_model.predict(sess, SEQ_LEN)
			np.savetxt(f'result.txt', np.array(fake_data).reshape((MAX_ATOMS, 2)), delimiter=';')

Run(sys.argv[1])
