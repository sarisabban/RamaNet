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
parser.add_argument('-db', '--DatasetBack', nargs='+', metavar='', help='Build the Backbone dataset')
parser.add_argument('-df', '--DatasetFrag', nargs='+', metavar='', help='Build the Fragment dataset')
parser.add_argument('-ds', '--DatasetSeq' , nargs='+', metavar='', help='Build the Sequence dataset')
parser.add_argument('-tb', '--TrainBack'  , action='store_true'  , help='Train the Backbone neural network')
parser.add_argument('-tf', '--TrainFrag'  , action='store_true'  , help='Train the Fragment neural network')
parser.add_argument('-ts', '--TrainSeq'   , action='store_true'  , help='Train the Sequence neural network')

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
				structure = Bio.PDB.PDBParser(QUIET=True).get_structure(TheName,
																		InFile)
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
				structure = Bio.PDB.PDBParser(QUIET=True).get_structure('X',
																	TheFile)
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
			structure = Bio.PDB.PDBParser(QUIET=True).get_structure('X',
																TheFile)
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
		HPCfile.write('''thefile=$(awk -v "line=${}" 'NR == line {}' ../PDB.list)\n'''.format('{PBS_ARRAY_INDEX}', '{ print; exit }'))
		HPCfile.write('{}/main/source/bin/relax.default.linuxgccrelease -relax:thorough -nstruct 100 -database {}/main/database -s $thefile'.format(path, path))
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
		with open('./dataset_PS.csv', 'w') as headPS:
			headPS.write(headerPS+'\n')
		# Setup dataset header for distance matrices
		headerCM = ['PDB_ID']
		for r in range(1, 150+1):
			for c in range(1, 150+1):
				headerCM.append(',aa{}_aa{}'.format(r, c))
		headerCM = ''.join(headerCM)
		with open('./dataset_CM.csv', 'w') as headCM:
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
				with open('dataset_PS.csv', 'a') as PSdata:
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
							if atom1-atom2 <= 12:
								CM.append(str(atom1-atom2))
							else:
								CM.append(str(0))
						except:
							CM.append(str(0))
				assert len(CM) == 22500
				ContactMap = ','.join(CM)
				TheCMLine = '{},{}\n'.format(File, ContactMap)
				with open('dataset_CM.csv', 'a') as CMdata:
					CMdata.write(TheCMLine)
			except: pass
	def VectorisePSCM(self, PS_file='dataset_PS.csv',
						CM_file='dataset_CM.csv',
						C_MAX=12,
						fp=np.float64):
		'''
		This function vectorises the datasets, normalises them, as well as
		constructs the final tensor and export the result as a serial.
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
		with h5py.File('PS+CM.hdf5', 'w') as data:
			dataset = data.create_dataset('default', data=dataset)
		# IMPORT WITH THIS COMMAND:
		#with h5py.File('PS+DM.hdf5', 'r') as data: dataset=data['default'][()]
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
						if p < 0: p = p + 360
						phi.append(p)
						s = pose.psi(r+1)
						if s < 0: s = s + 360
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
							assert	len(aa) == len(ss) == len(phi)\
							== len(psi) == len(sasa) == math.sqrt(len(ctmp))
							for AA,SS,P,S,SASA in zip(aa,ss,phi,psi,sasa):
								info.append('{},{},{},{},{}'\
								.format(AA, SS, P, S, SASA))
							Info = ','.join(info)
							with open('./AsPSa.csv', 'a') as data:
								data.write(File + ',' + Info + '\n')
							with open('lengths.txt', 'a') as length:
								length.write(str(len(aa))+'\n')
							for x in ctmp:
								m.append('{}'.format(x))
							M = ','.join(m)
							with open('./M.csv', 'a') as data:
								data.write(File + ',' + M + '\n')
							os.system('mv {} ./Completed'.format(TheFile))
						except: passos.system('mv {} ./Error_NotEqual'.format(TheFile))
				else: os.system('mv {} ./Error_Small'.format(TheFile))
			except: passos.system('mv {} ./Error_Broken'.format(TheFile))
	def Fill(self, filename):
		''' Fills missing .csv table spaces with zeros '''
		with open(filename) as f:
			with open(filename, 'a') as F:
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
	def HeaderAsPSaM(self, length=745, choice='AsPSa'):
		'''
		Constructs a .csv header and completes the dataset. To find the value of
		the largest structure run: sort -nk 1 lengths.txt
		'''
		header = ['PDB_ID']
		if choice == 'AsPSa':
			for i in range(1, length+1):
				header.append(',aa_{},ss_{},phi_{},psi_{},sasa_{}'\
				.format(i, i, i, i, i))
			header = ''.join(header)
			with open('./AsPSa.csv', 'r') as data:
				with open('./dataset_AsPSa.csv', 'w') as head:
					head.write(header+'\n')
					for line in data:
						head.write(line)
		elif choice == 'M':
			for r in range(1, length+1):
				for c in range(1, length+1):
					header.append(',aa{}_aa{}'.format(r, c))
			header = ''.join(header)
			with open('./M.csv', 'r') as data:
				with open('./dataset_M.csv', 'w') as head:
					head.write(header+'\n')
					for line in data:
						head.write(line)
	def build(self, switches='', directory='PDBDatabase'):
		if len(switches) == 15:
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
			if switch[14] == '1': self.DatasetPSCM('PDBCleaned')
			if switch[15] == '1': self.C_Max('dataset_CM.csv')
			if switch[16] == '1': self.VectorisePSCM()
		else: print('\x1b[31m[-] Error\x1b[33m: wrong string length\x1b[0m')

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

class Seq():
	def Database(self, TempDIR, FinalDIR):
		'''
		Downloads the entire PDB database from https://www.wwpdb.org/
		moves all files into one directory, then uncompresses all the files
		Generates a directory which contains all .PDB structure files
		'''
		os.system('rsync -rlpt -v -z --delete --port=33444 rsync.wwpdb.org::ftp/data/structures/divided/pdb/ ./{}'.format(TempDIR))
		os.mkdir(FinalDIR)
		filelist = os.listdir(TempDIR)
		print('\x1b[32m[+] Download complete\x1b[0m')
		print('\x1b[32m[+] Moving files\x1b[0m')
		for directories in tqdm.tqdm(filelist):
			files = os.listdir('{}/{}'.format(TempDIR, directories))
			for afile in files:
				location = ('{}/{}/{}'.format(TempDIR, directories, afile))
				os.rename(location, '{}/{}'.format(FinalDIR, afile))
		os.system('rm -r ./{}'.format(TempDIR))
	def Extract(self, directory):
		'''
		Extracts all the .ent.gz files and separate all chains and save them into
		seperate .pdb files. Replaces each .ent.gz file with the .pdb file of each
		chain
		'''
		current = os.getcwd()
		pdbfilelist = os.listdir(directory)
		io = Bio.PDB.PDBIO()
		os.chdir(directory)
		print('\x1b[32m[+] Extracting files\x1b[0m')
		for TheFile in tqdm.tqdm(pdbfilelist):
			try:
				TheName = TheFile.split('.')[0].split('pdb')[1].upper()
				InFile = gzip.open(TheFile, 'rt')
				structure = Bio.PDB.PDBParser(QUIET=True).get_structure(TheName, InFile)
				count = 0
				for chain in structure.get_chains():
					io.set_structure(chain)
					io.save(structure.get_id()+'_'+chain.get_id()+'.pdb')
				os.remove(TheFile)
			except Exception as TheError:
				print('\x1b[31m[-] Failed to extract\t{}\x1b[33m{}\x1b[0m'.format(TheFile.upper(), str(TheError)))
				os.remove(TheFile)
		os.chdir(current)
	def NonProtein(self, directory):
		''' Remove non-protein structures '''
		current = os.getcwd()
		pdbfilelist = os.listdir(directory)
		os.chdir(directory)
		print('\x1b[32m[+] Deleting none-protein structures\x1b[0m')
		for TheFile in tqdm.tqdm(pdbfilelist):
			structure = Bio.PDB.PDBParser(QUIET=True).get_structure('X', TheFile)
			ppb = Bio.PDB.Polypeptide.PPBuilder()
			Type = ppb.build_peptides(structure, aa_only=True)
			if Type == []:
				os.remove(TheFile)
			else:
				continue
		os.chdir(current)
	def Break(self, directory):
		''' Remove structures with a broken (non-continuous) chains '''
		current = os.getcwd()
		pdbfilelist = os.listdir(directory)
		os.chdir(directory)
		print('\x1b[32m[+] Removing structures with non-continuous chains\x1b[0m')
		for TheFile in tqdm.tqdm(pdbfilelist):
			structure = Bio.PDB.PDBParser(QUIET=True).get_structure('X', TheFile)
			ppb = Bio.PDB.Polypeptide.PPBuilder()
			Type = ppb.build_peptides(structure, aa_only=True)
			try:
				x = Type[1]
				os.remove(TheFile)
			except:
				continue
		os.chdir(current)
	def Renumber(self, directory):
		''' Renumber structures starting at 1 '''
		current = os.getcwd()
		pdbfilelist = os.listdir(directory)
		os.chdir(directory)
		print('\x1b[32m[+] Renumbering structures\x1b[0m')
		for TheFile in tqdm.tqdm(pdbfilelist):
			pdb = open(TheFile , 'r')
			PDB = open(TheFile + 'X' , 'w')
			count = 0
			num = 0
			AA2 = None
			for line in pdb:
				count += 1
				AA1 = line[23:27]
				if not AA1 == AA2:
					num += 1
				final_line = line[:7] + '{:4d}'.format(count) + line[11:17] + line[17:21] + 'A' + '{:4d}'.format(num) + line[26:]
				AA2 = AA1
				PDB.write(final_line)
			PDB.close()
			os.remove(TheFile)
			os.rename(TheFile + 'X' , TheFile)
		os.chdir(current)

def DatasetAsPSaM(directory):
	'''
	Compile a dataset of each residue's amino acid identify, secondary
	structure, phi angle, psi angle, solvent accessible surface area as
	a .csv file and the contact map as a separate .csv file. to be run
	after clean() on the ./cleaned directory and identifying the number
	of residuis of the largest structure
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
					if p < 0: p = p + 360
					phi.append(p)
					s = pose.psi(r+1)
					if s < 0: s = s + 360
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
						assert	len(aa) == len(ss) == len(phi)\
						== len(psi) == len(sasa) == math.sqrt(len(ctmp))
						for AA,SS,P,S,SASA in zip(aa,ss,phi,psi,sasa):
							info.append('{},{},{},{},{}'\
							.format(AA, SS, P, S, SASA))
						Info = ','.join(info)
						with open('./AsPSa.csv', 'a') as data:
							data.write(File + ',' + Info + '\n')
						with open('lengths.txt', 'a') as length:
							length.write(str(len(aa))+'\n')
						for x in ctmp:
							m.append('{}'.format(x))
						M = ','.join(m)
						with open('./M.csv', 'a') as data:
							data.write(File + ',' + M + '\n')
						os.system('mv {} ./Completed'.format(TheFile))
					except: os.system('mv {} ./Error_NotEqual'.format(TheFile))
			else: os.system('mv {} ./Error_Small'.format(TheFile))
		except: os.system('mv {} ./Error_Broken'.format(TheFile))
	def Fill(self, filename):
		''' Fills missing .csv table spaces with zeros '''
		with open(filename) as f:
			with open(filename, 'a') as F:
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
	def Header(self, length=745, choice='AsPSa'):
		'''
		Constructs a .csv header and completes the dataset. To find the value of
		the largest structure run: sort -nk 1 lengths.txt
		'''
		header = ['PDB_ID']
		if choice == 'AsPSa':
			for i in range(1, length+1):
				header.append(',aa_{},ss_{},phi_{},psi_{},sasa_{}'\
				.format(i, i, i, i, i))
			header = ''.join(header)
			with open('./AsPSa.csv', 'r') as data:
				with open('./dataset_AsPSa.csv', 'w') as head:
					head.write(header+'\n')
					for line in data:
						head.write(line)
		elif choice == 'M':
			for r in range(1, length+1):
				for c in range(1, length+1):
					header.append(',aa{}_aa{}'.format(r, c))
			header = ''.join(header)
			with open('./M.csv', 'r') as data:
				with open('./dataset_M.csv', 'w') as head:
					head.write(header+'\n')
					for line in data:
						head.write(line)

class BACKBONE():
	def gan(self, filename1=None, filename2=None, choice='generate'):
		'''
		A generative adversarial neural network that generates novel unnatural
		protein backbone topologies. This network uses the phi and psi angles
		as well as a distance matrix as protein structure features
		'''
		lrG = 0.001
		lrD = 0.001
		nodeG = 6
		nodeD = 5
		moment = 0.8
		alpha = 0.2
		drop = 0.25
		kernel = 3
		stride = 2
		latent = 128
		batchs = 64
		epochs = 1
		C_MAX = 12
		try:
			data = pd.read_csv(filename1)
			phi = data[data.columns[1::2]].values
			psi = data[data.columns[2::2]].values
			phi /= 180
			psi /= 180
			phi -= 1
			psi -= 1
			data = pd.read_csv(filename2)
			cm = data[data.columns[1:]].values
			cm = np.reshape(cm, (-1, 150, 150))
			cm /= (C_MAX/2)
			cm -= 1
			CM = cm
			PS = np.array([phi, psi])
			PS = np.swapaxes(PS, 0, 2)
			PS = np.swapaxes(PS, 0, 1)
			dataset = np.concatenate([PS, CM], axis=2)
			dataset = np.reshape(dataset,
			(-1, dataset.shape[1], dataset.shape[2], 1))
			PS, CM, cm, phi, psi = [], [], [], [], []
			sklearn.utils.shuffle(dataset)
			shape = dataset.shape[1:]
			print(dataset.shape)
		except: shape = (150, 152, 1)
		G = Sequential()
		G.add(Dense(2**(nodeG+1) * 75 * 38, activation='relu',input_dim=latent))
		G.add(Reshape((75, 38, 2**(nodeG+1))))
		G.add(UpSampling2D(size=(1, 2)))
		G.add(Conv2D(2**(nodeG+1), kernel_size=kernel, padding='same'))
		G.add(BatchNormalization(momentum=moment))
		G.add(Activation('relu'))
		G.add(UpSampling2D())
		G.add(Conv2D(2**(nodeG+0), kernel_size=kernel, padding='same'))
		G.add(BatchNormalization(momentum=moment))
		G.add(Activation('relu'))
		G.add(Conv2D(1, kernel_size=kernel, padding='same'))
		G.add(Activation('tanh'))
		D = Sequential()
		D.add(Conv2D(2**(nodeD+0), kernel_size=kernel, strides=stride, input_shape=shape, padding='same'))
		D.add(LeakyReLU(alpha=alpha))
		D.add(Dropout(drop))
		D.add(Conv2D(2**(nodeD+1), kernel_size=kernel, strides=stride, padding='same'))
		D.add(ZeroPadding2D(padding=((0, 1), (0, 1))))
		D.add(BatchNormalization(momentum=moment))
		D.add(LeakyReLU(alpha=alpha))
		D.add(Dropout(drop))
		D.add(Conv2D(2**(nodeD+2), kernel_size=kernel, strides=stride, padding='same'))
		D.add(BatchNormalization(momentum=moment))
		D.add(LeakyReLU(alpha=alpha))
		D.add(Dropout(drop))
		D.add(Conv2D(2**(nodeD+3), kernel_size=kernel, strides=stride-1, padding='same'))
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
		if choice == 'train':
			Epc, DTy, DFy, GNy = [], [], [], []
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
				if epoch % (k+1) == 0:
					g_loss = AM.train_on_batch(X_noise, y_true)
				GN_loss = round(float(g_loss[0]), 3)
				Epc.append(epoch)
				DTy.append(DT_loss)
				DFy.append(DF_loss)
				GNy.append(GN_loss)
				Verb =	'Epoch: {:6d} [DT {:.7f}][DF {:.7f}][G {:.7f}]'\
						.format(epoch, DT_loss, DF_loss, GN_loss)
				print(Verb)
			G.save_weights('weights.h5')
			return(Epc, DTy, DFy, GNy)
		if choice == 'generate':
			try: G.load_weights('weights.h5')
			except: print('Missing file: weights.h5'); exit()
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
			return(P, S, C)
	def train(self, filename=None):
		''' Train the neural network '''
		self.gan('train', filename)
	def generate(self, structures=1):
		''' Generate structures and fold them '''
		for i in range(1, structures+1):
			while True:
				P, S, C = self.gan()
				try: self.fold(P, S, C)
				except: continue
				if self.SQM('backbone.pdb')[1] == True:
					os.system('mv backbone.pdb {}.pdb'.format(str(i)))
				else: os.remove('backbone.pdb')
	def SQM(self, filename):
		'''
		Structure Quality Metric:
		Calculates the ratio of helices and sheets to loops, the radius of
		gyration, and the percent of amino acids comprising the structure core.
		Then averages their values. Returns a value between 0.0-1.0 where 1.0
		is a good structure.
		'''
		parser = Bio.PDB.PDBParser()
		structure = parser.get_structure('{}'.format(filename), filename)
		dssp = Bio.PDB.DSSP(structure[0], filename, acc_array='Wilke')
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
	def fold(self, P, S, C):
		''' Folds a structure using phi/psi angles and contact map '''
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
			print('[-] Generated structure not satisfactory')
			raise ValueError('Unsatisfactory structure')
		size = pose.residues.__len__()
		with open('constraints.cst', 'w') as thefile:
			for a in range(1, size+1):
				for A in range(1, size+1):
					if C[a][A] !=0:
						line = 'AtomPair CA {} CA {} GAUSSIANFUNC {} 1.0\n'\
						.format(a, A, C[a][A])
						thefile.write(line)
		con = pyrosetta.rosetta.protocols.constraint_movers.ConstraintSetMover()
		con.constraint_file('constraints.cst')
		con.add_constraints(True)
		con.apply(pose)
		scorefxn = get_fa_scorefxn()
		score_manager = pyrosetta.rosetta.core.scoring.ScoreTypeManager()
		atom_pair_constraint = score_manager.score_type_from_name('atom_pair_constraint')
		rama_prepro = score_manager.score_type_from_name('rama_prepro')
		scorefxn.set_weight(atom_pair_constraint, 5)
		scorefxn.set_weight(rama_prepro, 5)
		relax = pyrosetta.rosetta.protocols.relax.FastRelax()
		relax.set_scorefxn(scorefxn)
		os.remove('turnicated.pdb')
		os.remove('constraints.cst')
		relax.apply(pose)
		pose.dump_pdb('backbone.pdb')

class FRAGMENT():
	''' A neural network that generates 3-mer and 9-mer fragments'''
	def vectorise(self, filename='Fragments.csv', nx=1452):
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
	def fold(self, p, s, o):
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
	def fragments(self, aa='AAA', ss='HHH', p=[], s=[], o=[]):
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
	def picking(self, aa='MSSRSELLLEKF', ss='LLLHHHHHHHHH', size=3):
		''' Fragment picking '''
		mer3, mer9 = [], []
		for i in range(1, size+1):
			p, s, o = lstm(choice='predict', aa=aa, ss=ss)
			III, IX = self.fragments(aa, ss, p, s, o)
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
	def lstm(self, X='X.hdf5', Y='Y.hdf5', choice='train', aa=[], ss=[]):
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

class SEQUENCE():
	'''  '''
	pass

def main():
	if args.DatasetBack:
		DB = Dataset()
		DB.build(sys.argv[2])
	elif args.DatasetFrag:  #### ADD TO READ ME
		DF = Vall()
		DF.vall()
	elif args.DatasetSeq: #### ADD TO READ ME
		DS = Dataset()
		DS.Database('DATABASE', 'PDBDatabase')
		DS.Extract('PDBDatabase')
		DS.NonProtein('PDBDatabase')
		DS.Break('PDBDatabase')
		DS.Renumber('PDBDatabase')
		DS.DatasetAsPSaM('cln')
		DS.Header(280, 'AsPSa')
		DS.Fill('dataset_header.csv')
	elif args.TrainBack:
		print('\x1b[33m[.] Training...\x1b[0m')
		BB = BACKBONE()
		BB.train('dataset_PS.csv', 'dataset_DM.csv')
		print('\x1b[32m[+] Training done\x1b[0m')

	 #### ADD TO READ ME GENERATE BACKBONE

	elif args.TrainFrag:  #### ADD TO READ ME
		F = FRAGMENT()
		F.lstm()
	 #### ADD TO READ ME GENERATE Fragments
	elif args.TrainSeq:  #### ADD TO READ ME
		pass
	 #### ADD TO READ ME GENERATE sequence
if __name__ == '__main__': main()
