import os
import sys
import Bio.PDB
import numpy as np
import pandas as pd
import tensorflow as tf
from pyrosetta import *
from pyrosetta.toolbox import *
init()

def LSTM_GAN(choice):
	TRAIN_DATA_FILE = './PS_Helix_500.csv'	# Dataset location
	NUM_EPOCHS = 3000						# Number of training epochs
	MAX_ATOMS = 150							# Maximum protein chain length
	SEQ_LEN = MAX_ATOMS * 2					# Total prediction sequence length (2 angles per atom)
	N_STRUCTURES = 1						# Number of structures to predict
	def FoldPDB_PS(data):
		size = int(len(data[0]))
		Vs = list()
		for numb in range(size): Vs.append('V')
		sequence = ''.join(Vs)
		pose = pose_from_sequence(sequence)
		PHI = data[0]
		PSI = data[1]
		count = 1
		for P, S in zip(PHI, PSI):
			pose.set_phi(count, float(P))
			pose.set_psi(count, float(S))
			count += 1
		atom = 1
		scorefxn = get_fa_scorefxn()
		relax = pyrosetta.rosetta.protocols.relax.FastRelax(scorefxn)
		pose.dump_pdb('temp1.pdb')
		structure = Bio.PDB.PDBParser().get_structure('temp', 'temp.pdb')
		dssp = Bio.PDB.DSSP(structure[0], 'temp.pdb', acc_array='Wilke')
		ppb = Bio.PDB.Polypeptide.PPBuilder()
		chain = ppb.build_peptides(structure, aa_only=False)[0]
		SS = []
		for aa in dssp:
			if aa[2] == 'G' or aa[2] == 'H' or aa[2] == 'I': SSname = 'H'
			elif aa[2] == 'B' or aa[2] == 'E': SSname = 'S'
			else: SSname = 'L'
			SS.append(SSname)
		# Adjust End
		for i in enumerate(reversed(SS)):
			if i[1] != 'L':
				num = i[0]
				break
		for model in structure:
			for chain in model:
				for i in reversed(range(150-num, 150+1)):
					chain.detach_child((' ', i, ' '))
		io = Bio.PDB.PDBIO()
		io.set_structure(structure)
		io.save('temp2.pdb')
		os.remove('temp1.pdb')
		pose = pose_from_pdb('temp2.pdb')
		# Just relax once is enough
		relax.apply(pose)
		# Simulated annealing relax (not nessesary)
		'''
		pose_R = Pose()
		for i in range(20):
			pose_R.assign(pose)
			score_B = scorefxn(pose)
			relax.apply(pose_R)
			score_A = scorefxn(pose_R)
			if score_A < score_B:
				pose.assign(pose_R)
		'''
		pose.dump_pdb('Backbone.pdb')
		os.remove('temp2.pdb')
	def Filter(TheFile):
		'''
		A function that filters protein structures
		'''
		structure = Bio.PDB.PDBParser().get_structure('{}'.format(TheFile), TheFile)
		dssp = Bio.PDB.DSSP(structure[0], TheFile, acc_array='Wilke')
		ppb = Bio.PDB.Polypeptide.PPBuilder()
		chain = ppb.build_peptides(structure, aa_only=False)[0]
		choice = True
		SS = []
		CST = []
		SASA = []
		for aa in dssp:
			if aa[2] == 'G' or aa[2] == 'H' or aa[2] == 'I':
				SSname = 'H'
			elif aa[2] == 'B' or aa[2] == 'E':
				SSname = 'S'
			else:
				SSname = 'L'
			SS.append(SSname)
			if   aa[1]=='A' : sasa=129*(aa[3])
			elif aa[1]=='V' : sasa=174*(aa[3])
			elif aa[1]=='I' : sasa=197*(aa[3])
			elif aa[1]=='L' : sasa=201*(aa[3])
			elif aa[1]=='M' : sasa=224*(aa[3])
			elif aa[1]=='P' : sasa=159*(aa[3])
			elif aa[1]=='Y' : sasa=263*(aa[3])
			elif aa[1]=='F' : sasa=240*(aa[3])
			elif aa[1]=='W' : sasa=285*(aa[3])
			elif aa[1]=='R' : sasa=274*(aa[3])
			elif aa[1]=='N' : sasa=195*(aa[3])
			elif aa[1]=='C' : sasa=167*(aa[3])
			elif aa[1]=='Q' : sasa=225*(aa[3])
			elif aa[1]=='E' : sasa=223*(aa[3])
			elif aa[1]=='G' : sasa=104*(aa[3])
			elif aa[1]=='H' : sasa=224*(aa[3])
			elif aa[1]=='K' : sasa=236*(aa[3])
			elif aa[1]=='S' : sasa=155*(aa[3])
			elif aa[1]=='T' : sasa=172*(aa[3])
			elif aa[1]=='D' : sasa=193*(aa[3])
			if sasa <= 15 and (SSname == 'H' or SSname == 'S'):
				layer = 'C'
			elif 15 < sasa < 60 and (SSname == 'H' or SSname == 'S'):
				layer = 'B'
			elif sasa >= 60 and (SSname == 'H' or SSname == 'S'):
				layer = 'S'
			if sasa <= 25 and SSname == 'L':
				layer = 'C'
			elif 25 < sasa < 40 and SSname == 'L':
				layer = 'B'
			elif sasa >= 40 and SSname == 'L':
				layer = 'S'
			SASA.append(layer)
			residue1 = chain[0]
			residue2 = chain[aa[0]-1]
			atom1 = residue1['CA']
			atom2 = residue2['CA']
			CST.append(atom1-atom2)
		# Secondary structure filter
		Hs = [w.replace('L', '.') for w in SS]
		Hs = [w.replace('S', '.') for w in Hs]
		Hs = ''.join(Hs).split('.')
		Hs = [item for item in Hs if item]
		Hnum = len(Hs)
		Ss = [w.replace('L', '.') for w in SS]
		Ss = [w.replace('H', '.') for w in Ss]
		Ss = ''.join(Ss).split('.')
		Ss = [item for item in Ss if item]
		Snum = len(Ss)
		Ls = [w.replace('H', '.') for w in SS]
		Ls = [w.replace('S', '.') for w in Ls]
		Ls = ''.join(Ls).split('.')
		Ls = [item for item in Ls if item]
		Lnum = len(Ls)
		H = SS.count('H')
		S = SS.count('S')
		L = SS.count('L')
		if len(SS) < 80:
			choice = False
		if H+S < L:
			choice = False
		Surface = SASA.count('S')
		Boundary = SASA.count('B')
		Core = SASA.count('C')
		percent = (Core*100)/(Surface+Boundary+Core)
		if percent < 15:
			choice = False
		MaxCST = max(CST)
		if MaxCST > 88:
			choice = False
		return(choice)
	class ModelConfig(object):
		def __init__(self):
			self.num_layers = 1			# Number of LSTM layers
			self.rnn_size = 64			# Number of LSTM units
			self.hidden_size = 32		# Number of hidden layer units
			self.num_mixtures = 4
			self.batch_size = 4
			self.num_steps = 10
			self.dropout_rate = 0.5		# Dropout rate
			self.learning_rate = 0.001	# Learning rate
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
	def Run(trn_prd):
		train_config = ModelConfig()
		train_config.learning_rate = 0.0003
		test_config = ModelConfig()
		test_config.batch_size = 1
		test_config.num_steps = 1
		# Train
		if trn_prd == 'train':
			os.makedirs('weights')
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
				saver.save(sess, f'./weights/weights.ckpt')
			true_data = data[0]
		# Predict
		elif trn_prd == 'predict':
			ckpt_path = f'./weights/weights.ckpt'
			for structure in range(N_STRUCTURES):
				reset_session_and_model()
				with tf.Session() as sess:
					test_model = MDNModel(test_config, True)
					test_model.is_training = False
					sess.run(tf.global_variables_initializer())
					saver = tf.train.Saver()
					saver.restore(sess, ckpt_path)
					fake_data = test_model.predict(sess, SEQ_LEN)
				np.savetxt(f'prediction.txt', np.array(fake_data).reshape((MAX_ATOMS, 2)), delimiter=';')
	if choice == 'train':
		Run(choice)
	elif choice == 'predict':
		while True:
			Run(choice)
			newfile = open('prediction.txt', 'r')
			phiout = []
			psiout = []
			for line in newfile:
				line = line.strip().split(';')
				phiout.append(float(line[0]))
				psiout.append(float(line[1]))
			phiout = [x*360.0 for x in phiout]
			psiout = [x*360.0 for x in psiout]
			data = (phiout, psiout)
			FoldPDB_PS(data)
			os.remove('prediction.txt')
			if Filter('Backbone.pdb'):
				break
			else:
				os.remove('Backbone.pdb')

LSTM_GAN(sys.argv[1])
