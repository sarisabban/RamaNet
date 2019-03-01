import os
import sys
import keras
import numpy as np
import pandas as pd

def LSTMGAN_PS(choice, dataset):
	'''
	An LSTM Generative Adverserial Neural Network that will learn the
	structure of ideal proteins given their phi, psi angles
	(the dataPS.csv dataset). Then it generates novel angles and from
	random noise that will fold into a novel protein backbone.
	'''
	#Network values
	shape = (150, 2)
	latent = 100
	batchs = 4
	epochs = 3
	maxlen = 150
	seqlen = maxlen * 2
	layers = 1
	size = 64
	hidden = 32
	mixtures = 4
	steps = 10
	dropout = 0.5
	lr = 0.001
	# Import data
	data = pd.read_csv(dataset, ';')
	PHI = data[data.columns[2::3]].values
	PSI = data[data.columns[3::3]].values
	# Normalise
	PHI /= 360
	PSI /= 360
	# Make tensor (examples, residues, angles)
	DATA = np.array([PHI, PSI])
	DATA = np.swapaxes(DATA, 0, 2)
	DATA = np.swapaxes(DATA, 0, 1)
	# Generate sequences and next angles
	for i in range(0, maxlen, steps):
		seq = DATA[i:i+maxlen]
		ang = DATA[i+1:i+maxlen+1]
		sequences = np.concatenate(seq)
		nextangle = np.concatenate(ang)


	chars = '1'
	#Discriminator
	D = keras.models.Sequential()
	D.add(keras.layers.LSTM(layers, input_shape=(maxlen, len(chars))))
	D.add(keras.layers.Dense(len(chars), activation='softmax'))
	#Generator
	G = keras.models.Sequential()
	G.add(keras.layers.LSTM(layers, input_shape=(maxlen, len(chars))))
	G.add(keras.layers.Dense(len(chars), activation='softmax'))


	#Discriminator Model
	DM = keras.models.Sequential()
	DM.add(D)
	DM.compile(optimizer=keras.optimizers.Adam(lr), loss='binary_crossentropy', metrics=['accuracy'])
	#Adversarial Model
	AM = keras.models.Sequential()
	AM.add(G)
	AM.add(D)
	AM.compile(optimizer=keras.optimizers.Adam(lr), loss='binary_crossentropy', metrics=['accuracy'])
	'''
	if choice == 'train':
		#Training
		for epoch in range(epochs):
			#Generate a fake structures
			real = X[np.random.randint(0, X.shape[0], size=batchs)]
			noise = np.random.uniform(0.0, 1.0, size=[batchs, 100])
			fake = G.predict(noise)
			#Train discriminator
			x = np.concatenate((real, fake))
			y = np.ones([2*batchs, 1])
			y[batchs:, :] = 0
			d_loss = DM.train_on_batch(x, y)
			#Train adversarial
			y = np.ones([batchs, 1])
			a_loss = AM.train_on_batch(noise, y)
			D_loss = round(float(d_loss[0]), 3)
			D_accu = round(float(d_loss[1]), 3)
			A_loss = round(float(a_loss[0]), 3)
			print('{:7} [D loss: {:.3f}, accuracy: {:.3f}] [G loss: {:.3f}]'.format(epoch, D_loss, D_accu, A_loss))
			#Save Model
			G.save_weights('weights.h5')
	elif choice == 'generate':
		#Generate
		G.load_weights('weights.h5')
		noise = np.random.normal(0.5, 0.5, (1, 100))
		gen = G.predict(noise)
		gen = gen.reshape([300])
		gen = np.ndarray.tolist(gen)
		phiout = gen[0::2]		#[start:end:step]
		psiout = gen[1::2]		#[start:end:step]
		#Re-normalise
		phiout = [x*360.0 for x in phiout]
		psiout = [x*360.0 for x in psiout]
		return(phiout, psiout)


start_index = random.randint(0 , len(text) - maxlen - 1)	#Choose a random number from 0 to length of text - max charachter length - 1. This gives a number that is that start of a poisiton in the text
sentence = text[start_index : start_index + maxlen]			#Use that position to slice out a sentance from start_index position with a length of maxlen max charachter length. This gives a random sentance from the text
print(sentence)												#Print the stating sentance
for iter in range(400):										#Move 400 steps. Controls length of text
	#One-hot encode that sentance
	x_pred = numpy.zeros((1 , maxlen , len(chars)))			#Generate a tensor that has the same of (1 , max sentance length , number of available charachter), in our example: (1 , 40 , 58). The tensor is just filled with zeros, no other information
	for t , char in enumerate(sentence):					#Loop through the randomly generated sentance
		x_pred[0 , t , chars_indices[char]] = 1.0			#One-hot encode the randomly generated sentance (put a 1.0 for each charachter as available from the list of characters)
	#Use that tensor to make a prediction of the next charachter that comes after the randomly generated sentance
	preds = model.predict(x_pred , verbose = 0)[0]			#Returns a vector of shape (number of available charachter,) with the values of the probability of each charachter being the next charachter after the randomly generated sentance
	#Decode that character
	temperature = 0.2										#Temperature is used to make lower probabilities lower and higher probabilities higher (using Temperature < 1.0) or vise versa (using Temperature < 1.0). Calibrate until the best temperatures is achieved. Temperature of 1 does nothing
	preds = numpy.asarray(preds).astype('float64')			#Make sure all tensor values are float64
	preds = numpy.log(preds) / temperature					#Log each tensor value and then divide each value by the temperature
	exp_preds = numpy.exp(preds)							#Turn each tensor value into an exponant
	preds = exp_preds / numpy.sum(exp_preds)				#Re-Normalise (all values add up to 1) by dividing the exponant values by the sum of all the values
	probas = numpy.random.multinomial(1 , preds , 1)		#Randomly choose one index based on probability (most times it will choose the index with the highest probability, but sometimes it will randomly choose a slightly lower one)
	next_index = numpy.argmax(probas)						#Choose the largest value's number location in the vector, which will correspond to the identify of the charachter from the charachter list "indices_chars"
	next_char = indices_chars[next_index]					#Find the value's corresponding charachter
	sentence = sentence[1 : ] + next_char					#Add new charachter to sentance and remove 1 charachter from start of the sentence to maintain its length
	sys.stdout.write(next_char)								#Print the generated characters from the neural network prediction
	sys.stdout.flush()										#Flush terminal buffer, this and the previous line allows for the characters to be printer like a type writer (one at a time)
print('\n--------------------')

		'''
LSTMGAN_PS('train', 'PSC_Helix_5.csv')
