import os
import re
import sys
import h5py
import time
import glob
import math
import tqdm
import gzip
import keras
import sherpa
import Bio.PDB
import datetime
import warnings
import argparse
import numpy as np
import pandas as pd
import tensorflow as tf
from pyrosetta import *
from pyrosetta.toolbox import *
from keras.optimizers import Adam
from keras.models import Sequential, Model
from keras.losses import BinaryCrossentropy
from keras.layers.convolutional import Conv2D
from keras.layers import Activation, ZeroPadding2D
from keras.layers.advanced_activations import LeakyReLU
from keras.layers import Input, Dense, Reshape, Flatten
from keras.layers import UpSampling2D, BatchNormalization
from keras.layers import Dropout, GlobalMaxPooling2D, Conv2DTranspose
init('-out:level 0')

#dataset = np.random.randint(low=0, high=1, size=(2, 150, 152, 1))
with h5py.File('PS+CM.hdf5', 'r') as data: dataset=data['default'][()]
dataset = np.reshape(dataset, (-1, 150, 152, 1))
shape = dataset.shape[1:]
os.makedirs('./saved', exist_ok=True)

def fold(P, S, C):
	P = np.ndarray.tolist(P)
	S = np.ndarray.tolist(S)
	size = int(len(P))
	Vs = []
	for numb in range(size): Vs.append('A')
	sequence = ''.join(Vs)
	pose = pose_from_sequence(sequence)
	for count, (phi, psi) in enumerate(zip(P, S)):
		pose.set_phi(count+1, float(phi))
		pose.set_psi(count+1, float(psi))
	pose.dump_pdb('angles.pdb')
	structure = Bio.PDB.PDBParser().get_structure('angles', 'angles.pdb')
	dssp = Bio.PDB.DSSP(structure[0], 'angles.pdb', acc_array='Wilke')
	#dssp = Bio.PDB.DSSP(structure[0], filename, acc_array='Wilke', dssp='/home/x_sabbans0a/.conda/envs/RoMLenv/bin/mkdssp') FOR USE WITH CONDA IN IBEX SYSTEM OF KAUST     RoMLenv is the conda env
	ppb = Bio.PDB.Polypeptide.PPBuilder()
	chain = ppb.build_peptides(structure, aa_only=False)[0]
	SS = []
	for aa in dssp:
		if aa[2] == 'G' or aa[2] == 'H' or aa[2] == 'I': SSname = 'H'
		elif aa[2] == 'B' or aa[2] == 'E': SSname = 'S'
		else: SSname = 'L'
		SS.append(SSname)
	try:
		for i in enumerate(reversed(SS)):
			if i[1] != 'L':
				num = i[0]
				break
		for model in structure:
			for chain in model:
				for i in reversed(range(150-num+1, 150+1)):
					chain.detach_child((' ', i, ' '))
		io = Bio.PDB.PDBIO()
		io.set_structure(structure)
		io.save('turnicated.pdb')
		os.remove('angles.pdb')
		pose = pose_from_pdb('turnicated.pdb')
	except:
		os.remove('angles.pdb')
		raise ValueError('No Secondary Structures')
	size = pose.residues.__len__()
	with open('constraints.cst', 'w') as thefile:
		for a in range(1, size+1):
			for A in range(1, size+1):
				line = 'AtomPair CA {} CA {} GAUSSIANFUNC {} 1.0\n'\
				.format(a, A, C[a][A])
				thefile.write(line)
	con = pyrosetta.rosetta.protocols.constraint_movers.ConstraintSetMover()
	con.constraint_file('constraints.cst')
	con.add_constraints(True)
	con.apply(pose)
	scorefxn = get_fa_scorefxn()
	score_manager = pyrosetta.rosetta.core.scoring.ScoreTypeManager()
	constraint = score_manager.score_type_from_name('atom_pair_constraint')
	scorefxn.set_weight(constraint, 5)
	relax = pyrosetta.rosetta.protocols.relax.FastRelax()
	relax.set_scorefxn(scorefxn)
	os.remove('turnicated.pdb')
	os.remove('constraints.cst')
	relax.apply(pose)
	pose.dump_pdb('backbone.pdb')

def SQM(filename):
	parser = Bio.PDB.PDBParser()
	structure = parser.get_structure('{}'.format(filename), filename)
	dssp = Bio.PDB.DSSP(structure[0], filename, acc_array='Wilke')
	#dssp = Bio.PDB.DSSP(structure[0], filename, acc_array='Wilke', dssp='/home/x_sabbans0a/.conda/envs/RoMLenv/bin/mkdssp') FOR USE WITH CONDA IN IBEX SYSTEM OF KAUST     RoMLenv is the conda env
	ppb = Bio.PDB.Polypeptide.PPBuilder()
	chain = ppb.build_peptides(structure, aa_only=False)[0]
	AminoAcid = {	'A':129, 'P':159, 'N':195, 'H':224,
					'V':174, 'Y':263, 'C':167, 'K':236,
					'I':197, 'F':240, 'Q':225, 'S':155,
					'L':201, 'W':285, 'E':223, 'T':172,
					'M':224, 'R':274, 'G':104, 'D':193}
	sec_struct = []
	SASA = []
	Ca_distances = []
	for aa in dssp:
		if   aa[2] == 'G' or aa[2] == 'H' or aa[2] == 'I': ss = 'H'
		elif aa[2] == 'B' or aa[2] == 'E':                 ss = 'S'
		elif aa[2] == 'S' or aa[2] == 'T' or aa[2] == '-': ss = 'L'
		sec_struct.append(ss)
		sasa = AminoAcid[aa[1]]*aa[3]
		if sasa <= 25:      sasa = 'C'
		elif 25 < sasa < 40:sasa = 'B'
		elif sasa >= 40:    sasa = 'S'
		SASA.append(sasa)
		residue1 = chain[0]
		residue2 = chain[aa[0]-1]
		atom1 = residue1['CA']
		atom2 = residue2['CA']
		Ca_distances.append(atom1-atom2)
	''' Secondary structure measurement '''
	H = len([x for x in sec_struct if x == 'H'])
	S = len([x for x in sec_struct if x == 'S'])
	L = len([x for x in sec_struct if x == 'L'])
	total = len(sec_struct)
	SS = (H+S)/total
	''' SASA measurement '''
	surface = len([x for x in SASA if x == 'S'])
	boundery = len([x for x in SASA if x == 'B'])
	core = len([x for x in SASA if x == 'C'])
	total = len(SASA)
	percent = (core*100)/total
	ratio = percent/20 #Cutoff point 20%
	if ratio > 1: Core = 1.0
	else: Core = ratio
	''' Radius of gyration measurement '''
	coord = list()
	mass = list()
	Structure = open(filename, 'r')
	for line in Structure:
		try:
			line = line.split()
			x = float(line[6])
			y = float(line[7])
			z = float(line[8])
			coord.append([x, y, z])
			if   line[-1] == 'C': mass.append(12.0107)
			elif line[-1] == 'O': mass.append(15.9994)
			elif line[-1] == 'N': mass.append(14.0067)
			elif line[-1] == 'S': mass.append(32.065)
		except: pass
	xm = [(m*i, m*j, m*k) for (i, j, k), m in zip(coord, mass)]
	tmass = sum(mass)
	rr = sum(mi*i + mj*j + mk*k for (i, j, k),(mi, mj, mk) in zip(coord,xm))
	mm = sum((sum(i)/tmass)**2 for i in zip(*xm))
	rg = math.sqrt(rr/tmass-mm)
	ratio = 15/rg #Cutoff point 15 angstroms
	if ratio > 1: Rg = 1.0
	else: Rg = ratio
	''' The metric '''
	Items = [SS, Rg, Core]
	TheSum = sum(Items)
	TheTotal = len(Items)
	TheMetric = TheSum/TheTotal
	''' The choice '''
	choice = True
	if len(sec_struct) < 80:   choice = False
	if H+S < L:                choice = False
	if percent < 15:           choice = False
	if max(Ca_distances) > 89: choice = False
	return((round(TheMetric, 5), choice))

def TheModel(lrG   = 0.001,
			lrD    = 0.001,
			nodeG  = 2,
			nodeD  = 2,
			moment = 0.99,
			alpha  = 0.2,
			drop   = 0.2,
			kernel = 2,
			stride = 2,
			latent = 100,
			batchs = 32,
			epochs = 10,
			C_MAX  = 12):
	G = Sequential()
	G.add(Dense(2**(nodeG+1) * 75 * 38, activation='relu',input_dim=latent))
	G.add(Reshape((75, 38, 2**(nodeG+1))))
	G.add(UpSampling2D(size=(1, 2)))
	G.add(Conv2D(2**(nodeG+1), kernel_size=(kernel, kernel), padding='same'))
	G.add(BatchNormalization(momentum=moment))
	G.add(Activation('relu'))
	G.add(UpSampling2D())
	G.add(Conv2D(2**(nodeG+0), kernel_size=(kernel, kernel), padding='same'))
	G.add(BatchNormalization(momentum=moment))
	G.add(Activation('relu'))
	G.add(Conv2D(1, kernel_size=(kernel, kernel), padding='same'))
	G.add(Activation('tanh'))
	D = Sequential()
	D.add(Conv2D(2**(nodeD+0), kernel_size=(kernel, kernel), strides=(stride, stride), input_shape=shape, padding='same'))
	D.add(LeakyReLU(alpha=alpha))
	D.add(Dropout(drop))
	D.add(Conv2D(2**(nodeD+1), kernel_size=(kernel, kernel), strides=(stride, stride), padding='same'))
	D.add(ZeroPadding2D(padding=((0, 1), (0, 1))))
	D.add(BatchNormalization(momentum=moment))
	D.add(LeakyReLU(alpha=alpha))
	D.add(Dropout(drop))
	D.add(Conv2D(2**(nodeD+2), kernel_size=(kernel, kernel), strides=(stride, stride), padding='same'))
	D.add(BatchNormalization(momentum=moment))
	D.add(LeakyReLU(alpha=alpha))
	D.add(Dropout(drop))
	D.add(Conv2D(2**(nodeD+3), kernel_size=(kernel, kernel), strides=(stride-1, stride-1), padding='same'))
	D.add(BatchNormalization(momentum=moment))
	D.add(LeakyReLU(alpha=alpha))
	D.add(Dropout(drop))
	D.add(Flatten())
	D.add(Dense(1, activation='sigmoid'))
	D.compile(optimizer=keras.optimizers.Adam(lrD), loss='binary_crossentropy', metrics=['accuracy'])
	z = keras.layers.Input(shape=(latent,))
	gen = G(z)
	D.trainable = False
	validity = D(gen)
	AM = keras.models.Model(z, validity)
	AM.compile(optimizer=keras.optimizers.Adam(lrG), loss='binary_crossentropy', metrics=['accuracy'])
	Epc = []
	DTy = []
	DFy = []
	GNy = []
	y_true = np.ones([batchs, 1])
	y_false = np.zeros([batchs, 1])
	k = 3
	for epoch in range(1, epochs+1):
		X_real = dataset[np.random.randint(0,
		dataset.shape[0],
		size=batchs)]
		X_noise = np.random.normal(0.0, 1.0, size=[batchs, latent])
		X_fake = G.predict(X_noise)
		dT_loss = D.train_on_batch(X_real, y_true)
		dF_loss = D.train_on_batch(X_fake, y_false)
		DT_loss = round(float(dT_loss[0]), 3)
		DF_loss = round(float(dF_loss[0]), 3)
		try: g_loss = [GNy[-1]]
		except: g_loss = [0]
		if epoch % (k+1) == 0: g_loss = AM.train_on_batch(X_noise, y_true)
		GN_loss = round(float(g_loss[0]), 3)
	metric = []
	for i in range(1, 10+1):
		noise = np.random.normal(0.0, 1.0, size=[1, latent])
		gen = G.predict(noise)
		P = gen[:,:,0]
		S = gen[:,:,1]
		C = gen[:,:,2:]
		P += 1
		S += 1
		C += 1
		P *= 180
		S *= 180
		C *= (C_MAX/2)
		P = np.reshape(P, (150,))
		S = np.reshape(S, (150,))
		C = np.reshape(C, (150, 150))
		try:
			fold(P, S, C)
			metric.append(SQM('backbone.pdb')[0])
			os.remove('backbone.pdb')
		except:
			metric.append(0)
	return(sum(metric)/len(metric))

parameters = [
	sherpa.Continuous('lrG',    [1e-2, 1e-5], 'log'),
	sherpa.Continuous('lrD',    [1e-2, 1e-5], 'log'),
	sherpa.Discrete(  'nodeG',  [4, 8]),
	sherpa.Discrete(  'nodeD',  [4, 5]),
	sherpa.Continuous('moment', [0.6, 0.99]),
	sherpa.Continuous('alpha',  [0.1, 0.5]),
	sherpa.Continuous('drop',   [0.0, 0.5]),
	sherpa.Choice(    'kernel', [2, 3, 5]),
	sherpa.Choice(    'stride', [2]),
	sherpa.Ordinal(   'latent', [128, 512]),
	sherpa.Ordinal(   'batchs', [32, 512]),
	sherpa.Choice(    'epochs', [1000, 3000, 5000, 7000, 9000, 10000])]

algorithm = sherpa.algorithms.PopulationBasedTraining(
	population_size=2,
	num_generations=5)

study = sherpa.Study(
	parameters=parameters,
	algorithm=algorithm,
	lower_is_better=False,
	disable_dashboard=True)

for trial in study:
	keras.backend.clear_session()
	generation = trial.parameters['generation']
	load_from  = trial.parameters['load_from']
	print('-'*55)
	print('Generation {} - Population {}'\
	.format(generation, trial.parameters['save_to']))
	metric = TheModel(	lrG    = trial.parameters['lrG'],
						lrD    = trial.parameters['lrD'],
						nodeG  = trial.parameters['nodeG'],
						nodeD  = trial.parameters['nodeD'],
						#moment = trial.parameters['moment'], # dtype here changes between 32 and 64 crashing the script therefore commented out
						alpha  = trial.parameters['alpha'],
						drop   = trial.parameters['drop'],
						kernel = trial.parameters['kernel'],
						stride = trial.parameters['stride'],
						latent = trial.parameters['latent'],
						batchs = trial.parameters['batchs'],
						epochs = trial.parameters['epochs'])
	print('Average Metric: {}'.format(metric))
	study.add_observation(	trial=trial,
							iteration=generation,
							objective=metric)
	study.finalize(trial=trial)
	study.save(output_dir='./saved')

print('\n=====BEST RESULTS=====')
results = study.get_best_result()
for key in results: print('{:>10}: {}'.format(key, results[key]))
