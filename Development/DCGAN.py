import keras
import numpy as np
import pandas as pd
from pyrosetta import *
from pyrosetta.toolbox import *
init()

def FoldPDB_PSC(data):
	size = int(len(data[0]))
	Vs = list()
	for numb in range(size):
		Vs.append('V')
	sequence = ''.join(Vs)
	pose = pose_from_sequence(sequence)
	PHI = data[0]
	PSI = data[1]
	CST = data[2]
	count = 1
	for P, S in zip(PHI, PSI):
		pose.set_phi(count, float(P))
		pose.set_psi(count, float(S))
		count += 1
	atom = 1
	for cst in CST:
		line = 'AtomPair CA 1 CA '+str(atom)+' GAUSSIANFUNC '+str(cst)+' 1.0\n'
		thefile = open('constraints.cst', 'a')
		thefile.write(line)
		thefile.close()
		atom += 1
	constraints = pyrosetta.rosetta.protocols.constraint_movers.ConstraintSetMover()
	constraints.constraint_file('constraints.cst')
	constraints.add_constraints(True)
	constraints.apply(pose)
	scorefxnCST = ScoreFunction()
	scorefxnCST.set_weight(pyrosetta.rosetta.core.scoring.ScoreType.atom_pair_constraint, 1.0)
	relaxCST = pyrosetta.rosetta.protocols.relax.FastRelax()
	relaxCST.set_scorefxn(scorefxnCST)
	relaxCST.constrain_relax_to_start_coords(True)
	relaxCST.constrain_coords(True)
	scorefxn = get_fa_scorefxn()
	relaxC = pyrosetta.rosetta.protocols.relax.FastRelax()
	relaxC.set_scorefxn(scorefxn)
	relaxC.constrain_relax_to_start_coords(True)
	relaxC.constrain_coords(True)
	relax = pyrosetta.rosetta.protocols.relax.FastRelax()
	relax.set_scorefxn(scorefxn)
#	relaxC.apply(pose)
	pose.dump_pdb('Backbone.pdb')
	os.remove('constraints.cst')

def DCGAN_PSC(choice, filename, CSTmax):
	#Network values
	shape = (150, 3)
	latent = 100
	batchs = 32
	epochs = 1000
	Disclr = 0.000001
	Advrlr = 0.00001
	# Import data
	data = pd.read_csv(filename, ';')
	# Convert data into numpy arrays
	phi = data[data.columns[2::3]].values
	psi = data[data.columns[3::3]].values
	cst = data[data.columns[4::3]].values
	# MinMax scaling
	phi /= 360
	psi /= 360
	cst /= float(CSTmax)
	# Make the tensor - shape (examples, residues, 3 channels 3 P S C)
	X = np.array([phi, psi, cst])	# Shape = (3, 82900, 150)
	X = np.swapaxes(X, 0, 2)		# Change shape to (150, 82900, 3)
	X = np.swapaxes(X, 0, 1)		# Change shape to (82900, 150, 3)
	#Discriminator - VGG19
	D = keras.models.Sequential()
	D.add(keras.layers.Conv1D(64, kernel_size=3, input_shape=shape))
	D.add(keras.layers.LeakyReLU(alpha=0.2))
	D.add(keras.layers.Conv1D(64, kernel_size=3))
	D.add(keras.layers.LeakyReLU(alpha=0.2))
	D.add(keras.layers.MaxPooling1D())
	D.add(keras.layers.Conv1D(128, kernel_size=3))
	D.add(keras.layers.LeakyReLU(alpha=0.2))
	D.add(keras.layers.Conv1D(128, kernel_size=3))
	D.add(keras.layers.LeakyReLU(alpha=0.2))
	D.add(keras.layers.MaxPooling1D())
	D.add(keras.layers.Conv1D(256, kernel_size=3))
	D.add(keras.layers.LeakyReLU(alpha=0.2))
	D.add(keras.layers.Conv1D(256, kernel_size=3))
	D.add(keras.layers.LeakyReLU(alpha=0.2))
	D.add(keras.layers.Conv1D(256, kernel_size=3))
	D.add(keras.layers.LeakyReLU(alpha=0.2))
	D.add(keras.layers.Conv1D(256, kernel_size=3))
	D.add(keras.layers.LeakyReLU(alpha=0.2))
#	D.add(keras.layers.MaxPooling1D())
	D.add(keras.layers.Conv1D(512, kernel_size=3))
	D.add(keras.layers.LeakyReLU(alpha=0.2))
	D.add(keras.layers.Conv1D(512, kernel_size=3))
	D.add(keras.layers.LeakyReLU(alpha=0.2))
	D.add(keras.layers.Conv1D(512, kernel_size=3))
	D.add(keras.layers.LeakyReLU(alpha=0.2))
	D.add(keras.layers.Conv1D(512, kernel_size=3))
	D.add(keras.layers.LeakyReLU(alpha=0.2))
#	D.add(keras.layers.MaxPooling1D())
	D.add(keras.layers.Conv1D(512, kernel_size=3))
	D.add(keras.layers.LeakyReLU(alpha=0.2))
	D.add(keras.layers.Conv1D(512, kernel_size=3))
	D.add(keras.layers.LeakyReLU(alpha=0.2))
	D.add(keras.layers.Conv1D(512, kernel_size=3))
	D.add(keras.layers.LeakyReLU(alpha=0.2))
	D.add(keras.layers.Conv1D(512, kernel_size=3))
	D.add(keras.layers.LeakyReLU(alpha=0.2))
	D.add(keras.layers.Flatten())
	D.add(keras.layers.Dense(1, activation='sigmoid'))
	#D.summary()
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
	#G.summary()
	#Discriminator Model
	DM = keras.models.Sequential()
	DM.add(D)
	DM.compile(optimizer=keras.optimizers.Adam(Disclr), loss='mean_squared_error', metrics=['accuracy'])
	#Adversarial Model
	AM = keras.models.Sequential()
	AM.add(G)
	AM.add(D)
	AM.compile(optimizer=keras.optimizers.Adam(Advrlr), loss='mean_squared_error', metrics=['accuracy'])
	if choice == 'train':
		#Training
		for epoch in range(epochs):
			#Generate a fake structure
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
			#if epoch % 100==0: print('{:7} [D loss: {:.7f}, accuracy: {:.7f}] [A loss: {:.7f}]'.format(epoch, D_loss, D_accu, A_loss))
			if epoch % 100==0: print('Epoch {:5} Discriminator loss: {:.5f} ||| Adversarial loss: {:.5f}'.format(epoch, D_loss, A_loss))
			#Save Model
			G.save_weights('weights.h5')
	elif choice == 'generate':
		#Generate
		G.load_weights('weights.h5')
		noise = np.random.normal(0.5, 0.5, (1, 100))
		gen = G.predict(noise)
		gen = gen.reshape([450])
		gen = np.ndarray.tolist(gen)
		phiout = gen[0::3]	#[start:end:step]
		psiout = gen[1::3]	#[start:end:step]
		cstout = gen[2::3]	#[start:end:step]
		#Re-normalise
		phiout = [x*360.0 for x in phiout]
		psiout = [x*360.0 for x in psiout]
		cstout = [x*float(CSTmax) for x in cstout]
		return(phiout, psiout, cstout)

DCGAN_PSC('train', 'PSC_Helix_5.csv', 88.731)
data = DCGAN_PSC('generate', 'PSC_Helix_500.csv', 88.731)
FoldPDB_PSC(data)
