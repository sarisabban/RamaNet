import os
import h5py
import keras
import random
import sklearn
import warnings
import numpy as np
from keras.optimizers import Adam
from keras.models import Sequential
from keras.preprocessing import sequence
from keras.layers import TimeDistributed, RepeatVector
from keras.layers import LSTM, Dense, Masking, Bidirectional
from pyrosetta import *
from pyrosetta.toolbox import *

def warn(*args, **kwargs): pass
warnings.warn = warn
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'
init('-out:level 0')
print('\x1b[36m--------------------------------------------------------\x1b[0m')

def Vall(filename='vall.jul19.2011', m=16800, nx=1490):
	'''
	Compile the PDB IDs, chains, phi, psi, omega, and SASA of all the structures
	from the Rosetta vall.jul19.2011 database into a .csv file
	'''
	assert os.path.isfile('./{}'.format(filename)),\
	'Make sure the vall.jul19.2011 file is in the same directory as this script'
	with open(filename, 'r') as f:
		with open('Fragments.csv', 'w') as F:
			header = ['PDBID,Chain']
			for i in range(1, nx+1):
				header.append(',AA_{},SS_{},P_{},S_{},O_{},SASA_{}'\
				.format(i, i, i, i, i, i))
			header = ''.join(header)
			F.write(header + '\n')
			for i in range(30): next(f)
			ID  = []
			CH  = []
			AA  = []
			SS  = []
			P   = []
			S   = []
			O   = []
			SASA= []
			ID_seen = set()
			for line in f:
				line = line.strip().split()
				if line[0] not in ID_seen:
					exp = []
					for aa, ss, p, s, o, sasa in zip(AA, SS, P, S, O, SASA):
						exp.append('{},{},{},{},{},{}'\
						.format(aa, ss, p, s, o, sasa))
					exp = ','.join(exp)
					if exp == '': pass
					else: F.write(ID + ',' + CH + ',' + exp + '\n')
					ID   = None
					CH   = None
					AA   = []
					SS   = []
					P    = []
					S    = []
					O    = []
					SASA = []
					ID_seen.add(line[0])
					ID = line[0][:4].upper()
					CH = line[0][-1].upper()
					AA.append(line[1])
					SS.append(line[2])
					P.append(line[14])
					S.append(line[15])
					O.append(line[16])
					SASA.append(line[19])
				else:
					ID = line[0][:4].upper()
					CH = line[0][-1].upper()
					AA.append(line[1])
					SS.append(line[2])
					P.append(line[14])
					S.append(line[15])
					O.append(line[16])
					SASA.append(line[19])
			exp = []
			for aa, ss, p, s, o, sasa in zip(AA, SS, P, S, O, SASA):
				exp.append('{},{},{},{},{},{}'\
				.format(aa, ss, p, s, o, sasa))
			exp = ','.join(exp)
			F.write(ID + ',' + CH + ',' + exp)

def vectorise(filename='Fragments.csv', nx=1452):
	''' Vectorises the dataset, normalises it, then serialises it '''
	# 1. Import data
	rows = len(open(filename).readlines()) - 1
	# 2. Generate a list of random number of rows
	lines = list(range(1, rows + 1))
	random.shuffle(lines)
	# 3. Open CSV file
	with open(filename, 'r') as File: all_lines_variable = File.readlines()
	PDBID, CHAIN, X, Y = [], [], [], []
	for i in lines:
		# 4. Import data line by line
		line = all_lines_variable[i]
		line = line.strip().split(',')
		if line[0] == '1OFD': continue # Causes an error
		aa   = np.array(line[2::6])
		ss   = np.array(line[3::6])
		p    = np.array(line[4::6])
		s    = np.array(line[5::6])
		o    = np.array(line[6::6])
		sasa = np.array(line[7::6])
		p    = np.array([float(i) for i in p])
		s    = np.array([float(i) for i in s])
		o    = np.array([float(i) for i in o])
		sasa = np.array([float(i) for i in sasa])
		# 5. Re-format data
		aa[aa=='A']     = 0
		aa[aa=='C']     = 1
		aa[aa=='D']     = 2
		aa[aa=='E']     = 3
		aa[aa=='F']     = 4
		aa[aa=='G']     = 5
		aa[aa=='H']     = 6
		aa[aa=='I']     = 7
		aa[aa=='K']     = 8
		aa[aa=='L']     = 9
		aa[aa=='M']     = 10
		aa[aa=='N']     = 11
		aa[aa=='P']     = 12
		aa[aa=='Q']     = 13
		aa[aa=='R']     = 14
		aa[aa=='S']     = 15
		aa[aa=='T']     = 16
		aa[aa=='V']     = 17
		aa[aa=='W']     = 18
		aa[aa=='Y']     = 19
		ss[ss=='L']     = 0
		ss[ss=='H']     = 1
		ss[ss=='E']     = 2
		p[p<0] = p[p<0] + 360
		s[s<0] = s[s<0] + 360
		o[o<0] = o[o<0] + 360
		aa = aa.astype(int)
		ss = ss.astype(int)
		# 6. Padding categories
		gap = nx - aa.size
		for pad in range(gap):
			aa = np.append(aa, -1)
			ss = np.append(ss, -1)
		# 7. One-hot encode amino acid sequences and secondary structures
		Aminos = []
		for x in aa:
			letter = [0 for _ in range(20)]
			if x != -1: letter[x] = 1
			Aminos.append(letter)
		Struct = []
		for x in ss:
			letter = [0 for _ in range(3)]
			if x != -1: letter[x] = 1
			Struct.append(letter)
		aa = np.array(Aminos)
		ss = np.array(Struct)
		# 8. Normalise data [min/max]
		p    = (p-0)/(360-0)
		s    = (s-0)/(360-0)
		o    = (o-0)/(360-0)
		sasa = (sasa-0)/(277-0)
		# 9. Padding values
		for pad in range(gap):
			p    = np.append(p,    0)
			s    = np.append(s,    0)
			o    = np.append(o,    0)
			sasa = np.append(sasa, 0)
		# 10. Expand axis
		p    = np.expand_dims(p,    axis=1)
		s    = np.expand_dims(s,    axis=1)
		o    = np.expand_dims(o,    axis=1)
		sasa = np.expand_dims(sasa, axis=1)
		# 11. Export
		featur = np.concatenate((aa, ss), axis=1)
		angles = np.concatenate((p, s, o), axis=1)
		PDBID.append(line[0])
		CHAIN.append(line[1])
		X.append(featur)
		Y.append(angles)
	PDBID = np.array(PDBID)
	CHAIN = np.array(CHAIN)
	PDBID = np.expand_dims(PDBID, axis=1)
	CHAIN = np.expand_dims(CHAIN, axis=1)
	X = np.array(X)
	Y = np.array(Y)
	print('X =', X.shape)
	print('Y =', Y.shape)
	# 12. Serialise tensors
	with h5py.File('Y.hdf5', 'w') as y:
		dset = y.create_dataset('default', data=Y)
	with h5py.File('X.hdf5', 'w') as x:
		dset = x.create_dataset('default', data=X)

def fold(p, s, o):
	''' Use the angle output of the LSTM network to fold a structure '''
	size = int(len(p))
	Vs = []
	for numb in range(size): Vs.append('V')
	sequence = ''.join(Vs)
	pose = pose_from_sequence(sequence)
	count = 1
	for P, S, O in zip(p, s, o):
		pose.set_phi(  count, float(P))
		pose.set_psi(  count, float(S))
		pose.set_omega(count, float(O))
		count += 1
	pose.dump_pdb('backbone.pdb')

def fragments(aa='AAA', ss='HHH', p=[], s=[], o=[]):
	''' Generate 3-mer and 9-mer fragments '''
	# 3-mer
	mer3 = []
	Pchunks3  = []
	Schunks3  = []
	Ochunks3  = []
	AAchunks3 = []
	SSchunks3 = []
	for i in range(int(len(aa)-2)):
		Pchunks3.append(p[i:i+3])
		Schunks3.append(s[i:i+3])
		Ochunks3.append(o[i:i+3])
		AAchunks3.append(aa[i:i+3])
		SSchunks3.append(ss[i:i+3])
	for item in zip(AAchunks3, SSchunks3, Pchunks3, Schunks3, Ochunks3):
		AA, SS, P, S, O = item[0], item[1], item[2], item[3], item[4]
		for AAs, SSs, Ps, Ss, Os in zip(AA, SS, P, S, O):
			Ps, Ss, Os = round(Ps, 3), round(Ss, 3), round(Os, 3)
			line = ' xxxx A    00 V {} {:8} {:8} {:8}'\
			.format(SSs, Ps, Ss, Os)
			mer3.append(line)
	# 9-mer
	mer9 = []
	Pchunks9  = []
	Schunks9  = []
	Ochunks9  = []
	AAchunks9 = []
	SSchunks9 = []
	for i in range(int(len(aa)-8)):
		Pchunks9.append(p[i:i+9])
		Schunks9.append(s[i:i+9])
		Ochunks9.append(o[i:i+9])
		AAchunks9.append(aa[i:i+9])
		SSchunks9.append(ss[i:i+9])
	for item in zip(AAchunks9, SSchunks9, Pchunks9, Schunks9, Ochunks9):
		AA, SS, P, S, O = item[0], item[1], item[2], item[3], item[4]
		for AAs, SSs, Ps, Ss, Os in zip(AA, SS, P, S, O):
			Ps, Ss, Os = round(Ps, 3), round(Ss, 3), round(Os, 3)
			line = ' xxxx A    00 V {} {:8} {:8} {:8}'\
			.format(SSs, Ps, Ss, Os)
			mer9.append(line)
	return(mer3, mer9)

def picking(aa='MSSRSELLLEKF', ss='LLLHHHHHHHHH', size=3):
	''' Fragment picking '''
	mer3, mer9 = [], []
	for i in range(1, size+1):
		p, s, o = lstm(choice='predict', aa=aa, ss=ss)
		III, IX = fragments(aa, ss, p, s, o)
		mer3.append(III)
		mer9.append(IX)
	with open('frags.200.9mers', 'w') as f:
		for n, N in enumerate(range(0, int(len(mer9[0])/9)*9, 9)):
			f.write('\n position: {:12} neighbors: {:9}\n'.format(n+1, 3))
			for P in range(len(mer9)):
				f.write('\n')
				for i in mer9[P][N:N+9]:
					f.write(i+'\n')
	with open('frags.200.3mers', 'w') as f:
		for n, N in enumerate(range(0, int(len(mer3[0])/3)*3, 3)):
			f.write('\n position: {:12} neighbors: {:9}\n'.format(n+1, 3))
			for P in range(len(mer3)):
				f.write('\n')
				for i in mer3[P][N:N+3]:
					f.write(i+'\n')

def lstm(X='X.hdf5', Y='Y.hdf5', choice='train', aa=[], ss=[]):
	''' Encoder-decoder sequence-to-sequence LSTM neural network '''
	try:
		with h5py.File(X, 'r') as x: X = x['default'][()]
		with h5py.File(Y, 'r') as y: Y = y['default'][()]
		shape = X.shape[1:]
	except: shape = (1452, 23)
	node1 = 150
	node2 = 150
	lr    = 0.001
	model = Sequential()
	model.add(Bidirectional(LSTM(node1), input_shape=shape))
	model.add(RepeatVector(1452))
	model.add(Bidirectional(LSTM(node2, return_sequences=True)))
	model.add(TimeDistributed(Dense(3)))
	model.compile(optimizer=Adam(lr=lr), loss='mean_squared_error')
	if choice == 'train':
		model.fit(X, Y, batch_size=32, epochs=1, verbose=1)
		model.save_weights('weights.h5')
	elif choice == 'predict':
		model.load_weights('weights.h5')
		aa = np.array([x for x in aa])
		ss = np.array([x for x in ss])
		aa[aa=='A'] = 0
		aa[aa=='C'] = 1
		aa[aa=='D'] = 2
		aa[aa=='E'] = 3
		aa[aa=='F'] = 4
		aa[aa=='G'] = 5
		aa[aa=='H'] = 6
		aa[aa=='I'] = 7
		aa[aa=='K'] = 8
		aa[aa=='L'] = 9
		aa[aa=='M'] = 10
		aa[aa=='N'] = 11
		aa[aa=='P'] = 12
		aa[aa=='Q'] = 13
		aa[aa=='R'] = 14
		aa[aa=='S'] = 15
		aa[aa=='T'] = 16
		aa[aa=='V'] = 17
		aa[aa=='W'] = 18
		aa[aa=='Y'] = 19
		ss[ss=='L'] = 0
		ss[ss=='H'] = 1
		ss[ss=='E'] = 2
		aa = aa.astype(int)
		ss = ss.astype(int)
		Aminos = []
		for x in aa:
			letter = [0 for _ in range(20)]
			if x != -1: letter[x] = 1
			Aminos.append(letter)
		Struct = []
		for x in ss:
			letter = [0 for _ in range(3)]
			if x != -1: letter[x] = 1
			Struct.append(letter)
		aa = np.array(Aminos)
		ss = np.array(Struct)
		X = np.concatenate((aa, ss), axis=1)
		X = np.expand_dims(X, axis=0)
		prediction = model.predict(X)[0]
		p = prediction[:,0]
		s = prediction[:,1]
		o = prediction[:,2]
		p = p * 360
		s = s * 360
		o = o * 360
		p = p.tolist()
		s = s.tolist()
		o = o.tolist()
		p = p[:len(aa)]
		s = s[:len(aa)]
		o = o[:len(aa)]
		return(p, s, o)

aa = 'MSSRSELLLEKFAEKIGIGSISFNENRLCSFAIDEIYYISLSDANDEYMMIYGVCGKFPTDNSNFALEILNANLWFAENGGPYLCYEAGAQSLLLALRFPLDDATPEKLENEIEVVVKSMENLYLVLHNQGITLKIEEISS'
ss = 'LLLHHHHHHHHHHHHHLLLLLLLLLLLLLEEEELLLLEEEEELLLLLEEEEEEEEEELLLLLHHHHHHHHHHHHHHHHLLLLEEEEELLLLEEEEEEEEELLLLLHHHHHHHHHHHHHHHHHHHHHHHLLLLLLLEEEELL'
#vectorise()
#lstm(choice='train')
#p, s, o = lstm(choice='predict', aa=aa, ss=ss)
#fold(p, s, o)
#picking(aa, ss, 200)
