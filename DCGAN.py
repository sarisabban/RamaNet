#!/usr/bin/python3

import keras
import numpy as np
import pandas as pd

def DCGAN():
	# Import data
	data = pd.read_csv('dataPSC.csv', ';')
	# Convert data into numpy arrays
	phi = data[data.columns[2::3]].values
	psi = data[data.columns[3::3]].values
	cst = data[data.columns[4::3]].values
	# MinMax scaling
	phi /= 360
	psi /= 360
	cst /= 207.801
	# Make the tensor - shape (examples, residues, 3 channels 3 P S C)
	X = np.array([phi, psi, cst])	# Shape = (3, 5187, 150)
	X = np.swapaxes(X, 0, 2)	# Change shape to (150, 5187, 3)
	X = np.swapaxes(X, 0, 1)	# Change shape to (5187, 150, 3)
	#Network values
	shape = (150, 3)
	latent = 100
	batchs = 32
	epochs = 1
	#Discriminator
	D = keras.models.Sequential()
	D.add(keras.layers.Conv1D(32, kernel_size=3, input_shape=shape))
	D.add(keras.layers.LeakyReLU(alpha=0.2))
	D.add(keras.layers.Conv1D(64, kernel_size=3))
	D.add(keras.layers.LeakyReLU(alpha=0.2))
	D.add(keras.layers.Conv1D(128, kernel_size=3))
	D.add(keras.layers.LeakyReLU(alpha=0.2))
	D.add(keras.layers.Conv1D(256, kernel_size=3))
	D.add(keras.layers.LeakyReLU(alpha=0.2))
	D.add(keras.layers.Flatten())
	D.add(keras.layers.Dense(1, activation='sigmoid'))
	D.summary()
	#Generator
	G = keras.models.Sequential()
	G.add(keras.layers.Dense(79*3, activation='relu', input_dim=latent))
	G.add(keras.layers.Reshape((79, 3)))
	G.add(keras.layers.Conv1D(128, kernel_size=3))
	G.add(keras.layers.Activation('relu'))
	G.add(keras.layers.UpSampling1D())
	G.add(keras.layers.Conv1D(64, kernel_size=3))
	G.add(keras.layers.Activation('relu'))
	G.add(keras.layers.Conv1D(3, kernel_size=3))
	G.add(keras.layers.Activation('tanh'))
	G.summary()
	#Discriminator Model
	DM = keras.models.Sequential()
	DM.add(D)
	DM.compile(optimizer=keras.optimizers.Adam(0.001), loss='binary_crossentropy', metrics=['accuracy'])
	#Adversarial Model
	AM = keras.models.Sequential()
	AM.add(G)
	AM.add(D)
	AM.compile(optimizer=keras.optimizers.Adam(0.001), loss='binary_crossentropy', metrics=['accuracy'])
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
		G.save_weights('DCGAN.h5')

DCGAN()
