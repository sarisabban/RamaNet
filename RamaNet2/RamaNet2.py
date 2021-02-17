#!/usr/bin/ python3

print('''\x1b[32m
██████╗  █████╗ ███╗   ███╗ █████╗ ███╗   ██╗███████╗████████╗
██╔══██╗██╔══██╗████╗ ████║██╔══██╗████╗  ██║██╔════╝╚══██╔══╝
██████╔╝███████║██╔████╔██║███████║██╔██╗ ██║█████╗     ██║
██╔══██╗██╔══██║██║╚██╔╝██║██╔══██║██║╚██╗██║██╔══╝     ██║
██║  ██║██║  ██║██║ ╚═╝ ██║██║  ██║██║ ╚████║███████╗   ██║
╚═╝  ╚═╝╚═╝  ╚═╝╚═╝     ╚═╝╚═╝  ╚═╝╚═╝  ╚═══╝╚══════╝   ╚═╝\x1b[35m
╔╦╗┌─┐  ┌┐┌┌─┐┬  ┬┌─┐  ╔═╗┬─┐┌─┐┌┬┐┌─┐┬┌┐┌  ╔╦╗┌─┐┌─┐┬┌─┐┌┐┌
 ║║├┤   ││││ │└┐┌┘│ │  ╠═╝├┬┘│ │ │ ├┤ ││││   ║║├┤ └─┐││ ┬│││
═╩╝└─┘  ┘└┘└─┘ └┘ └─┘  ╩  ┴└─└─┘ ┴ └─┘┴┘└┘  ═╩╝└─┘└─┘┴└─┘┘└┘
\u001b[31mAuthors:       \x1b[33mSari Sabban and Mikhail Markovsky
\u001b[31mDate:          \x1b[33m31-May-2017
\u001b[31mCorrespondace: \x1b[33msari.sabban@gmail.com
\u001b[31mURL:           \x1b[33mhttps://sarisabban.github.io/RamaNet
\x1b[36m---------------------------------------------------------\x1b[0m''')

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
import random
import sklearn
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

# Silence Tensorflow, Keras, and initialise PyRosetta
def warn(*args, **kwargs): pass
warnings.warn = warn
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'
init('-out:level 0')
print('\x1b[36m--------------------------------------------------------\x1b[0m')

# Setup arguments
parser = argparse.ArgumentParser(description='De Novo Protein Design Neural Network')
parser.add_argument('-d',  '--dataset',   nargs='+', metavar='', help='Build the Backbone or Sequence datasets')
parser.add_argument('-f',  '--frag',      action='store_true',   help='Build the Fragment dataset')
parser.add_argument('-tb', '--TrainBack', action='store_true',   help='Train the Backbone neural network')
parser.add_argument('-tf', '--TrainFrag', action='store_true',   help='Train the Fragment neural network')
parser.add_argument('-ts', '--TrainSeq',  action='store_true',   help='Train the Sequence neural network')
args = parser.parse_args()

class Dataset():
	''' Build a machine learning dataset of protein structures '''
	def Database(self, TempDIR, FinalDIR):
		'''
		Downloads the entire PDB database from https://www.wwpdb.org/
		moves all files into one directory, then uncompresses all the files
		Generates a directory which contains all .PDB structure files
		'''
		print('\x1b[33m[.] Downloading PDB database...\x1b[0m')
		web = 'rsync.wwpdb.org::ftp/data/structures/divided/pdb/'
		os.system('rsync -rlpt -q -v -z --delete --port=33444 {} {}'
		.format(web, TempDIR))
		print('\x1b[32m[+] Download complete\x1b[0m')
		os.mkdir(FinalDIR)
		filelist = os.listdir(TempDIR)
		print('\x1b[33m[.] Moving files...\x1b[0m')
		for directories in tqdm.tqdm(filelist):
			files = os.listdir('{}/{}'.format(TempDIR, directories))
			for afile in files:
				location = ('{}/{}/{}'.format(TempDIR, directories, afile))
				os.rename(location, '{}/{}'.format(FinalDIR, afile))
		os.system('rm -r ./{}'.format(TempDIR))
		print('\x1b[32m[+] Moving complete\x1b[0m')
	def Extract(self, directory):
		'''
		Extracts all the .ent.gz files and separate all chains and save them
		into seperate .pdb files. Replaces each .ent.gz file with the .pdb
		file of each chain
		'''
		print('\x1b[33m[.] Extracting files...\x1b[0m')
		current = os.getcwd()
		pdbfilelist = os.listdir(directory)
		io = Bio.PDB.PDBIO()
		os.chdir(directory)
		for TheFile in tqdm.tqdm(pdbfilelist):
			try:
				TheName = TheFile.split('.')[0].split('pdb')[1].upper()
				InFile = gzip.open(TheFile, 'rt')
				structure = Bio.PDB.PDBParser(QUIET=True)\
				.get_structure(TheName, InFile)
				count = 0
				for chain in structure.get_chains():
					io.set_structure(chain)
					io.save(structure.get_id()+'_'+chain.get_id()+'.pdb')
				os.remove(TheFile)
			except Exception as TheError:
				print('\x1b[31m[-] Failed to extract\t{}\x1b[33m: {}\x1b[0m'
						.format(TheFile.upper(), str(TheError)))
				os.remove(TheFile)
		os.chdir(current)
	def NonProtein(self, directory):
		''' Remove non-protein structures '''
		print('\x1b[33m[.] Deleting none-protein structures...\x1b[0m')
		current = os.getcwd()
		pdbfilelist = os.listdir(directory)
		os.chdir(directory)
		for TheFile in tqdm.tqdm(pdbfilelist):
			try:
				structure = Bio.PDB.PDBParser(QUIET=True)\
				.get_structure('X', TheFile)
				ppb = Bio.PDB.Polypeptide.PPBuilder()
				Type = ppb.build_peptides(structure, aa_only=True)
				if Type == []: os.remove(TheFile)
				else: continue
			except: os.remove(TheFile)
		os.chdir(current)
	def Size(self, directory, Size_From, Size_To):
		''' Remove structures not within defined size '''
		print('\x1b[33m[.] Removing unwanted structure sizes...\x1b[0m')
		current = os.getcwd()
		pdbfilelist = os.listdir(directory)
		os.chdir(directory)
		for TheFile in tqdm.tqdm(pdbfilelist):
			try:
				parser = Bio.PDB.PDBParser()
				structure = parser.get_structure('X', TheFile)
				model = structure[0]
				dssp = Bio.PDB.DSSP(model, TheFile, acc_array='Wilke')
				for aa in dssp: length = aa[0]
				if length >= int(Size_To) or length <= int(Size_From):
					os.remove(TheFile)
			except: print('\x1b[31m[-] Error in finding protein size\x1b[0m')
		os.chdir(current)
	def Break(self, directory):
		''' Remove structures with a broken (non-continuous) chains '''
		print('\x1b[33m[.] Removing non-continuous structures...\x1b[0m')
		current = os.getcwd()
		pdbfilelist = os.listdir(directory)
		os.chdir(directory)
		for TheFile in tqdm.tqdm(pdbfilelist):
			structure = Bio.PDB.PDBParser(QUIET=True)\
			.get_structure('X', TheFile)
			ppb = Bio.PDB.Polypeptide.PPBuilder()
			Type = ppb.build_peptides(structure, aa_only=True)
			try:
				x = Type[1]
				os.remove(TheFile)
			except: continue
		os.chdir(current)
	def Loops(self, directory, LoopLength):
		'''
		Remove structures that have loops that are larger than a
		spesific length
		'''
		print('\x1b[33m[.] Removing structures with long loops...\x1b[0m')
		current = os.getcwd()
		pdbfilelist = os.listdir(directory)
		os.chdir(directory)
		for TheFile in tqdm.tqdm(pdbfilelist):
			try:
				parser = Bio.PDB.PDBParser()
				structure = parser.get_structure('X', TheFile)
				model = structure[0]
				dssp = Bio.PDB.DSSP(model, TheFile, acc_array='Wilke')
				SS = list()
				for res in dssp:
					ss = res[2]
					if ss == '-' or ss == 'T' or ss == 'S': SS.append('L')
					else: SS.append('.')
				loops = ''.join(SS).split('.')
				loops = [item for item in loops if item]
				LargeLoop = None
				for item in loops:
					if len(item) <= LoopLength: continue
					else: LargeLoop = 'LargeLoop'
				if LargeLoop == 'LargeLoop': os.remove(TheFile)
				else: continue
			except: os.remove(TheFile)
		os.chdir(current)
	def Renumber(self, directory):
		''' Renumber structures starting at 1 '''
		print('\x1b[33m[.] Renumbering structures...\x1b[0m')
		current = os.getcwd()
		pdbfilelist = os.listdir(directory)
		os.chdir(directory)
		for TheFile in tqdm.tqdm(pdbfilelist):
			pdb = open(TheFile, 'r')
			PDB = open(TheFile+'X', 'w')
			count = 0
			num = 0
			AA2 = None
			for line in pdb:
				count += 1
				AA1 = line[23:27]
				if not AA1 == AA2: num += 1
				final_line =line[:7]+'{:4d}'.format(count)+line[11:17]+\
							line[17:21]+'A'+'{:4d}'.format(num)+line[26:]
				AA2 = AA1
				PDB.write(final_line)
			PDB.close()
			os.remove(TheFile)
			os.rename(TheFile+'X', TheFile)
		os.chdir(current)
	def Rg(self, directory, RGcutoff):
		''' Remove structures that are below the Raduis of Gyration's value '''
		print('\x1b[33m[.] Removing structure low Rg values...\x1b[0m')
		current = os.getcwd()
		pdbfilelist = os.listdir(directory)
		os.chdir(directory)
		for TheFile in tqdm.tqdm(pdbfilelist):
			mass = list()
			Structure = open(TheFile, 'r')
			for line in Structure:
				line = line.split()
				if line[0] == 'TER' or line[0] == 'END': continue
				else:
					if line[-1] == 'C':   mass.append(12.0107)
					elif line[-1] == 'O': mass.append(15.9994)
					elif line[-1] == 'N': mass.append(14.0067)
					elif line[-1] == 'S': mass.append(32.0650)
					elif line[-1] == 'H': mass.append(1.00794)
					else: continue
			coord = list()
			p = Bio.PDB.PDBParser()
			structure = p.get_structure('X', TheFile)
			for model in structure:
				for chain in model:
					for residue in chain:
						for atom in residue: coord.append(atom.get_coord())
			xm = [(m*i, m*j, m*k) for (i, j, k), m in zip(coord, mass)]
			tmass = sum(mass)
			rr = sum(mi*i + mj*j + mk*k for (i, j, k), (mi, mj, mk)\
			in zip(coord, xm))
			mm = sum((sum(i)/tmass)**2 for i in zip( * xm))
			rg = math.sqrt(rr/tmass-mm)
			if rg <= RGcutoff: os.remove(TheFile)
			else: continue
		os.chdir(current)
	def Clean(self, directory):
		''' Clean each structure within a directory '''
		print('\x1b[33m[.] Cleaning structures...\x1b[0m')
		os.mkdir('PDBCleaned')
		current = os.getcwd()
		pdbfilelist = os.listdir(directory)
		os.chdir(directory)
		for TheFile in tqdm.tqdm(pdbfilelist):
			CurFile = open(TheFile, 'r')
			NewFile = open('Clean-{}'.format(TheFile), 'a')
			for line in CurFile:
				if line.split()[0] == 'ATOM': NewFile.write(line)
			CurFile.close()
			NewFile.close()
			os.system('mv Clean-{} ../PDBCleaned'.format(TheFile))
		os.chdir(current)
	def Path(self, directory, path):
		''' Generate a file with the path to each file '''
		print('\x1b[33m[.] Generating paths...\x1b[0m')
		current = os.getcwd()
		pdbfilelist = os.listdir(directory)
		os.chdir(directory)
		PathFile = open('PDB.list', 'a')
		for TheFile in tqdm.tqdm(pdbfilelist):
			line = '{}/PDBCleaned/{}\n'.format(path, TheFile)
			PathFile.write(line)
		os.system('mv PDB.list ../')
		os.chdir(current)
	def RelaxHPC(self, path, cores):
		'''
		Generate a PBS job scheduler to perform each structure
		relax on a HPC
		'''
		HPCfile = open('relax.pbs', 'w')
		HPCfile.write('#!/bin/bash\n')
		HPCfile.write('#PBS -N Relax\n')
		HPCfile.write('#PBS -q fat\n')
		HPCfile.write('#PBS -l select=1:ncpus=1\n')
		HPCfile.write('#PBS -j oe\n')
		HPCfile.write('#PBS -J 1-{}\n'.format(str(cores)))
		HPCfile.write('cd $PBS_O_WORKDIR\n')
		HPCfile.write('mkdir PDBRelaxed\n')
		HPCfile.write('cd PDBRelaxed\n')
		HPCfile.write('''thefile=$(awk -v "line=${PBS_ARRAY_INDEX}"''')
		HPCfile.write(''''NR == line { print; exit }' ../PDB.list)\n''')
		HPCfile.write('{}/main/source/bin/'.format(path))
		HPCfile.write('relax.default.linuxgccrelease')
		HPCfile.write('-relax:thorough -nstruct 100 -database ')
		HPCfile.write('{}/main/database -s $thefile'.format(path))
		print('\x1b[32m[+] Generated HPC job submission file\x1b[0m')
	def Relax(self, directory):
		''' Relax each structure in a directory on a local computer '''
		print('\x1b[33m[.] Relaxing structures...\x1b[0m')
		os.mkdir('PDBRelaxed')
		current = os.getcwd()
		pdbfilelist = os.listdir(directory)
		os.chdir(directory)
		for TheFile in tqdm.tqdm(pdbfilelist):
			for i in range(1, 101):
				scorefxn = get_fa_scorefxn()
				relax = pyrosetta.rosetta.protocols.relax.FastRelax()
				relax.set_scorefxn(scorefxn)
				pose = pose_from_pdb(TheFile)
				relax.apply(pose)
				pose.dump_pdb('Relaxed{}-{}'.format(i, TheFile))
				os.system('mv Relaxed{}-{} ../PDBRelaxed'.format(i, TheFile))
		os.chdir(current)
	def C_Max(self, filename):
		''' Find the maximum value of the Distance Map in a dataset '''
		max_in_line = []
		with open(filename, 'r') as f:
			next(f)
			for line in f:
				line = line.strip().split(',')[1:]
				line = [float(item) for item in line]
				max_in_line.append(max(line))
			maximum = max(max_in_line)
			print('\x1b[32m[+] Contact Map maximum value: {}\x1b[0m'\
			.format(maximum))
			return(maximum)
	def DatasetPSCM(self, directory):
		'''
		Compile a dataset of each residue's phi and psi angles and another
		dataset of the contact map for each structure. This dataset is padded
		with zeros.
		'''
		a = 'Compiling phi and psi angles dataset'
		b = 'as well as a distance matrix dataset'
		text = a+b
		print('\x1b[32m{}\x1b[0m'.format(text))
		# Setup dataset header for angles
		headerPS = ['PDB_ID']
		for i in range(1, 150+1):
			headerPS.append(',phi_{},psi_{}'.format(i, i))
		headerPS = ''.join(headerPS)
		with open('./PS.csv', 'w') as headPS:
			headPS.write(headerPS+'\n')
		# Setup dataset header for distance matrices
		headerCM = ['PDB_ID']
		for r in range(1, 150+1):
			for c in range(1, 150+1):
				headerCM.append(',{}{}'.format(r, c))
		headerCM = ''.join(headerCM)
		with open('./CM.csv', 'w') as headCM:
			headCM.write(headerCM+'\n')
		for File in tqdm.tqdm(os.listdir(directory)):
			TheFile = '{}/{}'.format(directory, File)
			try:
				# Compile angles
				pose = pose_from_pdb(TheFile)
				phi = []
				psi = []
				for aa in range(len(pose.residues)):
					try:
						p = pose.phi(aa+1)
						s = pose.psi(aa+1)
						if p < 0: p = p+360
						if s < 0: s = s+360
						phi.append(p)
						psi.append(s)
					except: pass
				angles = []
				for P, S in zip(phi, psi):
					angles.append(str(round(P, 5))+','+str(round(S, 5)))
				assert len(phi) == len(psi)
				Angles = ','.join(angles)
				if len(angles) >= 150: AngLine = Angles
				else:
					addition = 150-len(angles)
					zeros = []
					for adds in range(addition): zeros.append('0.0,0.0')
					Zeros = ','.join(zeros)
					AngLine = '{},{}'.format(Angles, Zeros)
				ThePSLine = '{},{}\n'.format(File, AngLine)
				with open('PS.csv', 'a') as PSdata:
					PSdata.write(ThePSLine)
				#Compile contact map (Ca-Ca contact <= 12 angstroms)
				BIO = Bio.PDB.PDBParser(QUIET=True)
				structure = BIO.get_structure('X', TheFile)
				ppb = Bio.PDB.Polypeptide.PPBuilder()
				Type = ppb.build_peptides(structure, aa_only=False)
				model = Type
				chain = model[0]
				CM = []
				for aa1 in range(0, 150):
					for aa2 in range(0, 150):
						try:
							residue1 = chain[aa1]
							residue2 = chain[aa2]
							atom1 = residue1['CA']
							atom2 = residue2['CA']
							if atom1-atom2 <= 12: CM.append(str(atom1-atom2))
							else: CM.append(str(0))
						except:
							CM.append(str(0))
				assert len(CM) == 22500
				ContactMap = ','.join(CM)
				TheCMLine = '{},{}\n'.format(File, ContactMap)
				with open('CM.csv', 'a') as CMdata:
					CMdata.write(TheCMLine)
			except: pass
	def VectorisePSCM(self, PS_file='PS.csv',
						CM_file='CM.csv',
						C_MAX=12,
						fp=np.float64):
		'''
		This function vectorises the backbone PS and CM datasets, normalises
		them, combines them, as well as constructs the final tensor and
		export the result as a serial.
		'''
		# 1. Import a single row of PS dataset
		with open(PS_file) as PSf:
			next(PSf)
			P, S = [], []
			for line in PSf:
				# 2. Isolate different angles
				line = line.strip().split(',')
				p = [float(item) for item in line[1::2]]
				s = [float(item) for item in line[2::2]]
				assert len(p) == len(s)
				P.append(np.array(p, dtype=fp))
				S.append(np.array(s, dtype=fp))
		with open(CM_file) as CMf:
			next(CMf)
			CM = []
			for line in CMf:
				# 3. Isolate different points
				line = [float(item) for item in line.strip().split(',')[1:]]
				cm = np.reshape(line, (150, 150))
				CM.append(np.array(cm, dtype=fp))
		# 4. Construct PS matrices
		P = np.array(P)
		S = np.array(S)
		# 5. Normalise PS angles (min/max) [-1, 1]
		P /= 180
		S /= 180
		P -= 1
		S -= 1
		PS = np.array([P, S])
		PS = np.swapaxes(PS, 0, 2)
		PS = np.swapaxes(PS, 0, 1)
		# 6. Construct CM matrices
		CM = np.array(CM)
		# 7. Normalise CM contact map (min/max) [-1, 1]
		CM /= (C_MAX/2)
		CM -= 1
		# 8. Construct final dataset matrix
		dataset = np.concatenate([PS, CM], axis=2)
		# 9. Suffle dataset
		sklearn.utils.shuffle(dataset)
		# 10. Serialise tensors
		with h5py.File('PS+CM.h5', 'w') as data:
			dataset = data.create_dataset('default', data=dataset)
	def DatasetAsPSaM(self, directory):
		'''
		Compile a dataset of each residue's amino acid identify, secondary
		structure, phi angle, psi angle, solvent accessible surface area as
		a .csv file and the contact map as a separate .csv file. to be run
		after clean() on the ./cleaned directory, also outputs a file
		identifying the sizes of structures, so the largest value can be used
		with HeaderAsPSaM()
		'''
		os.makedirs('./Completed', exist_ok=True)
		os.makedirs('./Error_NotEqual', exist_ok=True)
		os.makedirs('./Error_Broken', exist_ok=True)
		os.makedirs('./Error_Small', exist_ok=True)
		for File in tqdm.tqdm(os.listdir(directory)):
			try:
				TheFile = '{}/{}'.format(directory, File)
				pose = pose_from_pdb(TheFile)
				DSSP = pyrosetta.rosetta.protocols.moves.DsspMover()
				DSSP.apply(pose)
				sasa_calc = pyrosetta.rosetta.core.scoring.sasa.SasaCalc()
				sasa_calc.calculate(pose)
				size = pose.total_residue()
				aa   = []
				ss   = []
				phi  = []
				psi  = []
				sasa = []
				info = []
				ctmp = []
				m    = []
				surf = list(sasa_calc.get_residue_sasa())
				for r  in range(size):
					if pose.residue(r+1).is_protein():
						aa.append(pose.sequence(r+1, r+1))
						ss.append(pose.secstruct(r+1))
						p = pose.phi(r+1)
						if p < 0: p = p+360
						phi.append(p)
						s = pose.psi(r+1)
						if s < 0: s = s+360
						psi.append(s)
						sasa.append(surf[r])
				for r  in range(0, size):
					for R  in range(0, size):
						if	pose.residue(r+1).is_protein() and\
							pose.residue(R+1).is_protein():
							CAr = pose.residue(r+1).xyz('CA')
							CAR = pose.residue(R+1).xyz('CA')
							CAr_CAR_vector = CAR-CAr
							Cont = CAr_CAR_vector.norm()
							if Cont <= 12: ctmp.append(Cont)
							else: ctmp.append(0)
				if len(aa) >= 50:
						try:
							assert len(aa) == len(ss) == len(phi)\
							== len(psi) == len(sasa) == math.sqrt(len(ctmp))
							for AA,SS,P,S,SASA in zip(aa,ss,phi,psi,sasa):
								info.append('{},{},{},{},{}'\
								.format(AA, SS, P, S, SASA))
							Info = ','.join(info)
							with open('./AsPSa_noheader_nofill.csv', 'a') as data:
								data.write(File + ',' + Info + '\n')
							with open('lengths.txt', 'a') as length:
								length.write(str(len(aa))+',')
							for x in ctmp:
								m.append('{}'.format(x))
							M = ','.join(m)
							with open('./M_noheader_nofill.csv', 'a') as data:
								data.write(File + ',' + M + '\n')
							os.system('mv {} ./Completed'.format(TheFile))
						except:
								passos.system('mv {} ./Error_NotEqual'\
								.format(TheFile))
				else: os.system('mv {} ./Error_Small'.format(TheFile))
			except: passos.system('mv {} ./Error_Broken'.format(TheFile))
	def HeaderAsPSaM(self, choice='AsPSa'):
		'''
		Constructs a .csv header and completes the dataset. To find the value of
		the largest structure run: sort -nk 1 lengths.txt
		'''
		with open('lengths.txt', 'r') as L:
			length = int(max(L.readlines()[0].strip().split(',')))
		header = ['PDB_ID']
		if choice == 'AsPSa':
			for i in range(1, length+1):
				header.append(',aa_{},ss_{},phi_{},psi_{},sasa_{}'\
				.format(i, i, i, i, i))
			header = ''.join(header)
			with open('./AsPSa_noheader_nofill.csv', 'r') as data:
				with open('./AsPSa_nofill.csv', 'w') as head:
					head.write(header+'\n')
					for line in data:
						head.write(line)
			os.remove('./AsPSa_noheader_nofill.csv')
		elif choice == 'M':
			for r in range(1, length+1):
				for c in range(1, length+1):
					header.append(',{}{}'.format(r, c))
			header = ''.join(header)
			with open('./M_noheader_nofill.csv', 'r') as data:
				with open('./M_nofill.csv', 'w') as head:
					head.write(header+'\n')
					for line in data:
						head.write(line)
			os.remove('./M_noheader_nofill.csv')
	def Fill(self, filename):
		''' Fills missing .csv table spaces with zeros '''
		name = filename.split('_')[0]
		with open(filename) as f:
			with open(name+'.csv', 'a') as F:
				first_line = f.readline()
				F.write(first_line)
				size = len(first_line.strip().split(','))
				for line in f:
					line = line.strip().split(',')
					gap = size - len(line)
					for zero in range(gap):
						line.append('0')
					new_line = ','.join(line)
					F.write(new_line + '\n')
		os.remove(filename)
	def VectoriseAsPSaM(self, filenameA='AsPSa.csv', filenameM='M.csv'):
		'''
		This function vectorises the backbone PS and CM datasets, normalises
		them, combines them, as well as constructs the final tensor and
		export the result as a serial.
		'''
		pass
	def build(self, switches='', directory='PDBDatabase'):
		if len(switches) == 20:
			switch = list(switches)
			if switch[0]  == '1': self.Database('DATABASE', directory)
			if switch[1]  == '1': self.Extract(directory)
			if switch[2]  == '1': self.NonProtein(directory)
			if switch[3]  == '1': self.Size(directory, 80, 150)
			if switch[4]  == '1': self.Break(directory)
			if switch[5]  == '1': self.Loops(directory, 10)
			if switch[6]  == '1': self.Renumber(directory)
			if switch[7]  == '1': self.Rg(directory, 15)
			########## --- HUMAN EYE FILTERING --- ##########
			if switch[8]  == '1': self.Clean(directory)
			if switch[9]  == '1': self.Path('PDBCleaned', '{PATH}')
			if switch[10] == '1': self.RelaxHPC('~/Rosetta', 829)
			if switch[11] == '1': self.Relax('PDBCleaned')
			if switch[12] == '1': self.DatasetAsPSaM('PDBCleaned')
			if switch[13] == '1': self.HeaderAsPSaM('AsPSa')
			if switch[14] == '1':
				self.HeaderAsPSaM('M')
				os.remove('lengths.txt')
			if switch[15] == '1':
				self.Fill('AsPSa_nofill.csv')
				self.Fill('M_nofill.csv')
			if switch[16] == '1': self.DatasetPSCM('PDBCleaned')
			if switch[17] == '1': self.C_Max('dataset_CM.csv')
			if switch[18] == '1': self.VectorisePSCM()
			if switch[18] == '1': self.VectoriseAsPSaM()
		else: print('\x1b[31m[-] Error\x1b[33m: Wrong string length\x1b[0m')

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

def Frag_vectorise(filename='Fragments.csv', nx=1452):
	''' Vectorises the fragments dataset, normalises it, then serialises it '''
	# 1. Import data
	rows = len(open(filename).readlines()) - 1
	# 2. Generate a list of random number of rows
	lines = list(range(1, rows + 1))
	random.shuffle(lines)
	# 3. Open CSV file
	with open(filename, 'r') as File: all_lines_variable = File.readlines()
	PDBID, CHAIN, X, Y = [], [], [], []
	for i in tqdm.tqdm(lines):
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
	with h5py.File('Frag_Y.h5', 'w') as y:
		dset = y.create_dataset('default', data=Y)
	with h5py.File('Frag_X.h5', 'w') as x:
		dset = x.create_dataset('default', data=X)

def SQM(filename):
	'''
	Structure Quality Metric:
	Calculates the ratio of helices and sheets to loops, the percent of amino
	acids comprising the structure core, and the radius of gyration as values
	between 0.0-1.0, it then averages the three values. Returns a value between
	0.0-1.0 where good structure >= 0.6
	'''
	parser = Bio.PDB.PDBParser()
	structure = parser.get_structure('{}'.format(filename), filename)
	dssp = Bio.PDB.DSSP(structure[0], filename, acc_array='Wilke')
	AminoAcid = {	'A':129, 'P':159, 'N':195, 'H':224,
					'V':174, 'Y':263, 'C':167, 'K':236,
					'I':197, 'F':240, 'Q':225, 'S':155,
					'L':201, 'W':285, 'E':223, 'T':172,
					'M':224, 'R':274, 'G':104, 'D':193}
	sec_struct = []
	SASA = []
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
	''' Secondary structure measurement '''
	H = len([x for x in sec_struct if x == 'H'])
	S = len([x for x in sec_struct if x == 'S'])
	L = len([x for x in sec_struct if x == 'L'])
	total = len(sec_struct)
	ratio = (H+S)/total
	limit = 1
	slope = 10
	bias  = 0.5
	SS = limit/(1+np.exp(slope*(bias-ratio)))
	''' SASA measurement '''
	surface = len([x for x in SASA if x == 'S'])
	boundery = len([x for x in SASA if x == 'B'])
	in_core = len([x for x in SASA if x == 'C'])
	total = len(SASA)
	percent = (in_core*100)/total
	Core = (2.50662/math.sqrt(2*(math.pi)))*math.exp(-((percent-30)**2)/100)
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
	rr = sum(mi*i + mj*j + mk*k for (i, j, k), (mi, mj, mk) in zip(coord, xm))
	mm = sum((sum(i)/tmass)**2 for i in zip(*xm))
	rg = math.sqrt(rr/tmass-mm)
	Rg = (2.50662/math.sqrt(2*(math.pi)))*math.exp(-((rg-12)**2)/40)
	''' The metric '''
	TheMetric = sum([SS, Core, Rg])/3
	return(round(TheMetric, 5))

class fold():
	''' Folds a protein structure given the phi/psi angles and contact map '''
	def __init__(self, Pa, Sa, CM):
		CM = np.reshape(CM, (150, 150))
		self.size = len([i for i in np.diag(CM, k=1) if i!=0])
		self.U = np.triu(CM, k=0)[:self.size, :self.size]
		self.L = np.tril(CM, k=0)[:self.size, :self.size]
		self.P = np.array(Pa)[:self.size]
		self.S = np.array(Sa)[:self.size]
	def upper_lower(self, side, name):
		''' Fold upper diagonal of contact map using the same phi/psi angles '''
		Vs = []
		for numb in range(self.size): Vs.append('A')
		sequence = ''.join(Vs)
		pose = pose_from_sequence(sequence)
		for count, (phi, psi) in enumerate(zip(self.P, self.S)):
			pose.set_phi(count+1, float(phi))
			pose.set_psi(count+1, float(psi))
		with open('constraints_{}.cst'.format(name), 'w') as thefile:
			for a in range(1, self.size):
				for A in range(1, self.size):
					if side[a][A] !=0:
						line = 'AtomPair CA {} CA {} GAUSSIANFUNC {} 1.0\n'\
						.format(a, A, side[a][A])
						thefile.write(line)
		con = pyrosetta.rosetta.protocols.constraint_movers.ConstraintSetMover()
		con.constraint_file('constraints_{}.cst'.format(name))
		con.add_constraints(True)
		con.apply(pose)
		scorefxn = get_fa_scorefxn()
		score_manager = pyrosetta.rosetta.core.scoring.ScoreTypeManager()
		atom_pair_constraint = score_manager.\
		score_type_from_name('atom_pair_constraint')
		rama_prepro = score_manager.score_type_from_name('rama_prepro')
		scorefxn.set_weight(atom_pair_constraint, 1)
		scorefxn.set_weight(rama_prepro, 1)
		relax = pyrosetta.rosetta.protocols.relax.FastRelax()
		relax.set_scorefxn(scorefxn)
		os.remove('constraints_{}.cst'.format(name))
		relax.apply(pose)
		pose.dump_pdb('backbone_{}.pdb'.format(name))
	def multi(self):
		''' Fold both upper and lower diagonals at the same time '''
		P1 = multiprocessing.Process(target=self.upper_lower,args=[self.U, 'U'])
		P2 = multiprocessing.Process(target=self.upper_lower,args=[self.L, 'L'])
		P1.start()
		P2.start()
		P1.join()
		P2.join()

################################################################################
############################### NEURAL NETWORKS ################################
################################################################################

class BACKBONE():
	''' A neural network that generates a protein backbone '''
	pass



class SEQUENCE():
	''' A neural network that generates a sequence for a protein backbone '''
	pass



class FRAGMENT():
	''' A neural network that generates 3-mer and 9-mer fragments '''
	pass




def main():
	if args.dataset:
		DB = Dataset()
		DB.build(sys.argv[2])
	elif args.frag:
		Vall()
		Frag_vectorise()
	elif args.TrainBack:
		print('\x1b[33m[.] Training...\x1b[0m')
		BB = BACKBONE()
		print('\x1b[32m[+] Training done\x1b[0m')
	elif args.TrainFrag:
		print('\x1b[33m[.] Training...\x1b[0m')
		FR = FRAGMENT()
		print('\x1b[32m[+] Training done\x1b[0m')
	elif args.TrainSeq:
		print('\x1b[33m[.] Training...\x1b[0m')
		SQ = SEQUENCE()
		print('\x1b[32m[+] Training done\x1b[0m')

if __name__ == '__main__': main()
