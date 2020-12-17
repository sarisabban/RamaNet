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
parser.add_argument('-d', '--dataset', nargs='+', metavar='', help='Build the dataset')
parser.add_argument('-t', '--train', action='store_true', help='Train the neural network')
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

























class GAN():
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

def main():
	if args.dataset:
		D = Dataset()
		D.build(sys.argv[2])
	elif args.train:
		print('\x1b[33m[.] Training...\x1b[0m')
		NN = GAN()
		NN.train('dataset_PS.csv', 'dataset_DM.csv')
		print('\x1b[32m[+] Training done\x1b[0m')

if __name__ == '__main__': main()
