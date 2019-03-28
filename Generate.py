#!/usr/bin/python

import os
import re
import bs4
import sys
import time
import glob
import keras
import Bio.PDB
import datetime
import requests
import argparse
import numpy as np
import pandas as pd
import urllib.request
import tensorflow as tf
from Bio import pairwise2
from pyrosetta import *
from pyrosetta.toolbox import *
init()

parser = argparse.ArgumentParser(description='De Novo Protein Design Neural Network')
parser.add_argument('-t', '--train', action='store_true', help='Train the neural network')
parser.add_argument('-f', '--fragments', action='store_true', help='Generate a structure and get its fragments from the Robetta server, you must specify a username')
args = parser.parse_args()

class RosettaDesign(object):
	def __init__(self, filename):
		''' Generate the resfile '''
		self.filename = filename
		parser = Bio.PDB.PDBParser()
		structure = parser.get_structure('{}'.format(filename), filename)
		dssp = Bio.PDB.DSSP(structure[0], filename, acc_array='Wilke')
		sasalist = []
		for x in dssp:
			if x[1] == 'A':
				sasa = 129*(x[3])
				if sasa <= 25:			sasa = 'C'
				elif 25 < sasa < 40:	sasa = 'B'
				elif sasa >= 40:		sasa = 'S'
			elif x[1] == 'V':
				sasa = 174*(x[3])
				if sasa <= 25:			sasa = 'C'
				elif 25 < sasa < 40:	sasa = 'B'
				elif sasa >= 40:		sasa = 'S'
			elif x[1] == 'I':
				sasa = 197*(x[3])
				if sasa <= 25:			sasa = 'C'
				elif 25 < sasa < 40:	sasa = 'B'
				elif sasa >= 40:		sasa = 'S'
			elif x[1] == 'L':
				sasa = 201*(x[3])
				if sasa <= 25:			sasa = 'C'
				elif 25 < sasa < 40:	sasa = 'B'
				elif sasa >= 40:		sasa = 'S'
			elif x[1] == 'M':
				sasa = 224*(x[3])
				if sasa <= 25:			sasa = 'C'
				elif 25 < sasa < 40:	sasa = 'B'
				elif sasa >= 40:		sasa = 'S'
			elif x[1] == 'P':
				sasa = 159*(x[3])
				if sasa <= 25:			sasa = 'C'
				elif 25 < sasa < 40:	sasa = 'B'
				elif sasa >= 40:		sasa = 'S'
			elif x[1] == 'Y':
				sasa = 263*(x[3])
				if sasa <= 25:			sasa = 'C'
				elif 25 < sasa < 40:	sasa = 'B'
				elif sasa >= 40:		sasa = 'S'
			elif x[1] == 'F':
				sasa = 240*(x[3])
				if sasa <= 25:			sasa = 'C'
				elif 25 < sasa < 40:	sasa = 'B'
				elif sasa >= 40:		sasa = 'S'
			elif x[1] == 'W':
				sasa = 285*(x[3])
				if sasa <= 25:			sasa = 'C'
				elif 25 < sasa < 40:	sasa = 'B'
				elif sasa >= 40:		sasa = 'S'
			elif x[1] == 'R':
				sasa = 274*(x[3])
				if sasa <= 25:			sasa = 'C'
				elif 25 < sasa < 40:	sasa = 'B'
				elif sasa >= 40:		sasa = 'S'
			elif x[1] == 'N':
				sasa = 195*(x[3])
				if sasa <= 25:			sasa = 'C'
				elif 25 < sasa < 40:	sasa = 'B'
				elif sasa >= 40:		sasa = 'S'
			elif x[1] == 'C':
				sasa = 167*(x[3])
				if sasa <= 25:			sasa = 'C'
				elif 25 < sasa < 40:	sasa = 'B'
				elif sasa >= 40:		sasa = 'S'
			elif x[1] == 'Q':
				sasa = 225*(x[3])
				if sasa <= 25:			sasa = 'C'
				elif 25 < sasa < 40:	sasa = 'B'
				elif sasa >= 40:		sasa = 'S'
			elif x[1] == 'E':
				sasa = 223*(x[3])
				if sasa <= 25:			sasa = 'C'
				elif 25 < sasa < 40:	sasa = 'B'
				elif sasa >= 40:		sasa = 'S'
			elif x[1] == 'G':
				sasa = 104*(x[3])
				if sasa <= 25:			sasa = 'C'
				elif 25 < sasa < 40:	sasa = 'B'
				elif sasa >= 40:		sasa = 'S'
			elif x[1] == 'H':
				sasa = 224*(x[3])
				if sasa <= 25:			sasa = 'C'
				elif 25 < sasa < 40:	sasa = 'B'
				elif sasa >= 40:		sasa = 'S'
			elif x[1] == 'K':
				sasa = 236*(x[3])
				if sasa <= 25:			sasa = 'C'
				elif 25 < sasa < 40:	sasa = 'B'
				elif sasa >= 40:		sasa = 'S'
			elif x[1] == 'S':
				sasa = 155*(x[3])
				if sasa <= 25:			sasa = 'C'
				elif 25 < sasa < 40:	sasa = 'B'
				elif sasa >= 40:		sasa = 'S'
			elif x[1] == 'T':
				sasa = 172*(x[3])
				if sasa <= 25:			sasa = 'C'
				elif 25 < sasa < 40:	sasa = 'B'
				elif sasa >= 40:		sasa = 'S'
			elif x[1] == 'D':
				sasa = 193*(x[3])
				if sasa <= 25:			sasa = 'C'
				elif 25 < sasa < 40:	sasa = 'B'
				elif sasa >= 40:		sasa = 'S'
			if x[2] == 'G' or x[2] == 'H' or x[2] == 'I':	ss = 'H'
			elif x[2] == 'B' or x[2] == 'E':				ss = 'S'
			elif x[2] == 'S' or x[2] == 'T' or x[2] == '-':	ss = 'L'
			sasalist.append((x[0], x[1], ss, sasa))
		resfile = open('resfile', 'a')
		resfile.write('NATRO\nSTART\n')
		for n, r, a, s in sasalist:
			if s == 'S' and a == 'L':
				line = '{} A PIKAA PGNQSTDERKH\n'.format(n)
				resfile.write(line)
			elif s == 'S' and a == 'H':
				line = '{} A PIKAA QEKH\n'.format(n)
				resfile.write(line)
			elif s == 'S' and a == 'S':
				line = '{} A PIKAA QTY\n'.format(n)
				resfile.write(line)
			elif s == 'B' and a == 'L':
				line = '{} A PIKAA AVILFYWGNQSTPDEKR\n'.format(n)
				resfile.write(line)
			elif s == 'B' and a == 'H':
				line = '{} A PIKAA AVILWQEKFM\n'.format(n)
				resfile.write(line)
			elif s == 'B' and a == 'S':
				line = '{} A PIKAA AVILFYWQTM\n'.format(n)
				resfile.write(line)
			elif s == 'C' and a == 'L':
				line = '{} A PIKAA AVILPFWM\n'.format(n)
				resfile.write(line)
			elif s == 'C' and a == 'H':
				line = '{} A PIKAA AVILFWM\n'.format(n)
				resfile.write(line)
			elif s == 'C' and a == 'S':
				line = '{} A PIKAA AVILFWM\n'.format(n)
				resfile.write(line)
		resfile.close()

	def __del__(self):
		''' Remove the resfile '''
		os.remove('resfile')
		try:	os.remove('fixbb.fasc')
		except:	os.remove('flxbb.fasc')

	def choose(self):
		''' Choose the lowest scoring structure '''
		try:	scorefile = open('fixbb.fasc', 'r')
		except:	scorefile = open('flxbb.fasc', 'r')
		score = 0
		name = None
		for line in scorefile:
			line = json.loads(line)
			score2 = line.get('total_score')
			if score2 < score:
				score = score2
				name = line.get('decoy')
		os.system('mv {} structure.pdb'.format(name))
		for f in glob.glob('f[il]xbb_*'): os.remove(f)

	def fixbb(self):
		'''
		Performs the RosettaDesign protocol to change a structure's
		amino acid sequence while maintaining a fixed backbone.
		Generates the structure.pdb file
		'''
		pose = pose_from_pdb(self.filename)
		starting_pose = Pose()
		starting_pose.assign(pose)
		scorefxn = get_fa_scorefxn()
		relax = pyrosetta.rosetta.protocols.relax.FastRelax()
		relax.set_scorefxn(scorefxn)
		packtask = standard_packer_task(pose)
		pyrosetta.rosetta.core.pack.task.parse_resfile(pose, packtask, 'resfile')
		fixbb = pyrosetta.rosetta.protocols.minimization_packing.PackRotamersMover(scorefxn, packtask, 10)
		job = PyJobDistributor('fixbb', 100, scorefxn)
		job.native_pose = starting_pose
		while not job.job_complete:
			pose.assign(starting_pose)
			relax.apply(pose)
			fixbb.apply(pose)
			relax.apply(pose)
			job.output_decoy(pose)
		self.choose()

	def flxbb(self):
		'''
		Performs the RosettaDesign protocol to change a structure's
		amino acid sequence while allowing for a flexible backbone.
		Generates the structure.pdb file
		'''
		pose = pose_from_pdb(self.filename)
		starting_pose = Pose()
		starting_pose.assign(pose)
		scorefxn = get_fa_scorefxn()
		relax = pyrosetta.rosetta.protocols.relax.FastRelax()
		relax.set_scorefxn(scorefxn)
		resfile = rosetta.core.pack.task.operation.ReadResfile('resfile')
		task = pyrosetta.rosetta.core.pack.task.TaskFactory()
		task.push_back(resfile)
		movemap = MoveMap()
		movemap.set_bb(True)
		movemap.set_chi(True)
		flxbb = pyrosetta.rosetta.protocols.denovo_design.movers.FastDesign()
		flxbb.set_task_factory(task)
		flxbb.set_movemap(movemap)
		flxbb.set_scorefxn(scorefxn)
		job = PyJobDistributor('flxbb', 100, scorefxn)
		job.native_pose = starting_pose
		while not job.job_complete:
			pose.assign(starting_pose)
			relax.apply(pose)
			flxbb.apply(pose)
			relax.apply(pose)
			job.output_decoy(pose)
		self.choose()

def Fragments(filename, username):
	'''
	Submits the pose to the Robetta server
	(http://www.robetta.org) for fragment generation that are
	used for the Abinitio folding simulation. Then measures the
	RMSD for each fragment at each position and chooses the
	lowest RMSD. Then averages out the lowest RMSDs. Then plots
	the lowest RMSD fragment for each positon.
	Generates the 3-mer file, the 9-mer file, the PsiPred file,
	the RMSD vs Position PDF plot with the averaged fragment
	RMSD printed in the plot
	'''
	#Make the 3-mer and 9-mer fragment files and the PSIPRED file using the Robetta server
	pose = pose_from_pdb(filename)
	sequence = pose.sequence()
	#Post
	web = requests.get('http://www.robetta.org/fragmentsubmit.jsp')
	payload = {	'UserName':		username,
			'Email':		'',
			'Notes':		'structure',
			'Sequence':		sequence,
			'Fasta':		'',
			'Code':			'',
			'ChemicalShifts':	'',
			'NoeConstraints':	'',
			'DipolarConstraints':	'',
			'type':			'submit'}
	session = requests.session()
	response = session.post('http://www.robetta.org/fragmentsubmit.jsp', data=payload , files=dict(foo='bar'))		
	for line in response:
		line = line.decode()
		if re.search('<a href="(fragmentqueue.jsp\?id=[0-9].*)">' , line):
			JobID = re.findall('<a href="(fragmentqueue.jsp\?id=[0-9].*)">' , line)
	JobURL = 'http://www.robetta.org/' + JobID[0]
	#Check
	ID = JobID[0].split('=')
	print('Job ID: ' + str(ID[1]))
	while True:
		Job = urllib.request.urlopen(JobURL)
		jobdata = bs4.BeautifulSoup(Job , 'lxml')
		status = jobdata.find('td', string='Status: ').find_next().text
		if status == 'Complete':
			print(datetime.datetime.now().strftime('%d %B %Y @ %H:%M') , 'Status:' , status)
			break
		else:
			print(datetime.datetime.now().strftime('%d %B %Y @ %H:%M') , 'Status:' , status)
			time.sleep(1800)
			continue
	#Download
	sequence = pose.sequence()
	fasta = open('structure.fasta' , 'w')
	fasta.write(sequence)
	fasta.close()
	time.sleep(1)
	os.system('wget http://www.robetta.org/downloads/fragments/' + str(ID[1])  + '/aat000_03_05.200_v1_3')
	os.system('wget http://www.robetta.org/downloads/fragments/' + str(ID[1])  + '/aat000_09_05.200_v1_3')
	os.system('wget http://www.robetta.org/downloads/fragments/' + str(ID[1])  + '/t000_.psipred_ss2')
	os.rename('aat000_03_05.200_v1_3' , 'frags.200.3mers')
	os.rename('aat000_09_05.200_v1_3' , 'frags.200.9mers')
	os.rename('t000_.psipred_ss2' , 'pre.psipred.ss2')
	#Calculate the best fragment's RMSD at each position
	frag = open('frags.200.9mers' , 'r')
	rmsd = open('temp.dat' , 'w')
	for line in frag:
		if line.lstrip().startswith('position:'):
			line = line.split()
			size = line[1]
	frag.close()
	count = 0
	for x in range (int(size)):
		count +=1
		#Get the pose and make a copy of it to apply changes to
		pose_copy = pyrosetta.Pose()
		pose_copy.assign(pose)
		#Setup frame list
		frames = pyrosetta.rosetta.core.fragment.FrameList()
		#Setup the 9-mer fragment (9-mer is better than 3-mer for this analysis)
		fragset = pyrosetta.rosetta.core.fragment.ConstantLengthFragSet(9)
		fragset.read_fragment_file('frags.200.9mers')
		fragset.frames(count , frames)
		#Setup the MoveMap
		movemap = MoveMap()
		movemap.set_bb(True)
		#Setup and apply the fragment inserting mover
		for frame in frames:
			for frag_num in range( 1 , frame.nr_frags() + 1 ):
				frame.apply(movemap , frag_num , pose_copy)
				#Measure the RMSD difference between the original pose and the new changed pose (the copy)
				RMSD = rosetta.core.scoring.CA_rmsd(pose , pose_copy)
				print(RMSD , '\t' , count)
				rmsd.write(str(RMSD) + '\t' + str(count) + '\n')
				#Reset the copy pose to original pose
				pose_copy.assign(pose)
	rmsd.close()
	#Analyse the RMSD file to get the lowest RMSD for each position
	data = open('RMSDvsPosition.dat' , 'w')
	lowest = {} 				#Mapping group number -> lowest value found
	for line in open('temp.dat'):
		parts = line.split()
		if len(parts) != 2:		#Only lines with two items on it
			continue
		first = float(parts[0])
		second = int(parts[1])
		if first == 0: 			#Skip line with 0.0 RMSD (this is an error from the 9-mer fragment file). I don't know why it happens
			continue
		if second not in lowest:
			lowest[second] = first
		else:
			if first < lowest[second]:
				lowest[second] = first
	for position, rmsd in lowest.items():
		#print(str(rmsd) + '\t' + str(position))
		data.write(str(position) + '\t' + str(rmsd) + '\n')
	data.close()
	#Calculate the average RMSD of the fragments
	data = open('RMSDvsPosition.dat' , 'r')
	value = 0
	for line in data:
		line = line.split()
		RMSD = float(line[1])
		value = value + RMSD
		count = int(line[0])
	Average_RMSD = round(value / count , 2)
	#Plot the results
	gnuplot = open('gnuplot_sets' , 'w')
	gnuplot.write("""
	reset\n
	set terminal postscript\n
	set output './plot_frag.pdf'\n
	set encoding iso_8859_1\n
	set term post eps enh color\n
	set xlabel 'Position'\n
	set ylabel 'RMSD (\\305)'\n
	set yrange [0:]\n
	set xrange [0:]\n
	set xtics auto\n
	set xtics rotate\n
	set grid front\n
	unset grid\n
	set title 'Fragment Quality'\n
	set key off\n
	set boxwidth 0.5\n
	set style fill solid\n
	set label 'Average RMSD = {}' at graph 0.01 , graph 0.95 tc lt 7 font 'curior 12'\n
	plot 'RMSDvsPosition.dat' with boxes\n
	exit
	""".format(str(Average_RMSD)))
	gnuplot.close()
	os.system('gnuplot < gnuplot_sets')
	os.remove('gnuplot_sets')
	os.remove('temp.dat')
	return(Average_RMSD)

def CSTMax(filename):
	'''
	find the minimum and maximum range of the constraints
	values of a dataset
	'''
	maxline = []
	data = open(filename, 'r')
	next(data)
	for line in data:
		line = line.strip().split(';')
		cst = []
		count = 1
		for item in line:
			if count < 450:
				count += 3
				cst.append(float(line[count]))
		maxline.append(max(cst))
	maximum = max(maxline)
	return(maximum)

def FoldPDB_PS(data):
	'''
	Fold a primary structure using the phi and psi torsion
	angles Generates the Backbone.pdb file
	'''
	#Generate a pose
	size = int(len(data[0]))
	Vs = list()
	for numb in range(size):
		Vs.append('V')
	sequence = ''.join(Vs)
	pose = pose_from_sequence(sequence)
	#Isolate each angle and constraint
	PHI = data[0]
	PSI = data[1]
	count = 1
	#Move amino acid angles
	for P, S in zip(PHI, PSI):
		pose.set_phi(count, float(P))
		pose.set_psi(count, float(S))
		count += 1
	atom = 1
	#Run FastRelax
	scorefxn = get_fa_scorefxn()
	relax = pyrosetta.rosetta.protocols.relax.FastRelax(scorefxn)
	relax.apply(pose)
	pose.dump_pdb('Backbone.pdb')

def FoldPDB_PSC(data):
	'''
	Fold a primary structure using the phi and psi torsion
	angles as well as the CA atom constraints. Generates
	the Backbone.pdb file
	'''
	#Generate a pose
	size = int(len(data[0]))
	Vs = list()
	for numb in range(size):
		Vs.append('V')
	sequence = ''.join(Vs)
	pose = pose_from_sequence(sequence)
	#Isolate each angle and constraint
	PHI = data[0]
	PSI = data[1]
	CST = data[2]
	count = 1
	#Move amino acid angles
	for P, S in zip(PHI, PSI):
		pose.set_phi(count, float(P))
		pose.set_psi(count, float(S))
		count += 1
	atom = 1
	#Write constraints file
	for cst in CST:
		line = 'AtomPair CA 1 CA '+str(atom)+' GAUSSIANFUNC '+str(cst)+' 1.0\n'
		thefile = open('constraints.cst', 'a')
		thefile.write(line)
		thefile.close()
		atom += 1
	#Add constraints option to pose
	constraints = pyrosetta.rosetta.protocols.constraint_movers.ConstraintSetMover()
	constraints.constraint_file('constraints.cst')
	constraints.add_constraints(True)
	constraints.apply(pose)
	#Setup score function with weight on only atom_pair_constraint
	scorefxnCST = ScoreFunction()
	scorefxnCST.set_weight(pyrosetta.rosetta.core.scoring.ScoreType.atom_pair_constraint, 1.0)
	#Setup constraint relax to bring atoms together
	relaxCST = pyrosetta.rosetta.protocols.relax.FastRelax()
	relaxCST.set_scorefxn(scorefxnCST)
	relaxCST.constrain_relax_to_start_coords(True)
	relaxCST.constrain_coords(True)
	#Setup normal FastRelax with constraints
	scorefxn = get_fa_scorefxn()
	relaxC = pyrosetta.rosetta.protocols.relax.FastRelax()
	relaxC.set_scorefxn(scorefxn)
	relaxC.constrain_relax_to_start_coords(True)
	relaxC.constrain_coords(True)
	#Setup normal FastRelax without constraints
	relax = pyrosetta.rosetta.protocols.relax.FastRelax()
	relax.set_scorefxn(scorefxn)
	#Run Relaxations
#	relaxCST.apply(pose)	#only brings structure together
	relaxC.apply(pose)		#Best on its own
#	relax.apply(pose)		#Best on its own
	pose.dump_pdb('Backbone.pdb')
	os.remove('constraints.cst')

def LSTM_GAN(choice):
	'''
	A neural network that designs a helical protein topology using phi/psi angels
	The neural network was chosen and optamised by Mikhail Markovsky.
	-----------------------------------------------------------------------------------
	This neural network architecture was inspired from https://github.com/nesl/sensegen 
	that has the following licence:

	All rights reserved Networked and Embedded Systems Lab (NESL), UCLA.
	Permission is hereby granted, free of charge, to any person obtaining a copy 
	of this software and associated documentation files (the "Software"), to deal 
	in the Software without restriction, including without limitation the rights 
	to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies 
	of the Software, and to permit persons to whom the Software is furnished to do 
	so, subject to the following conditions:
	The above copyright notice and this permission notice shall be included in all 
	copies or substantial portions of the Software.
	THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR 
	IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS 
	FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR 
	COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER 
	IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION 
	WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
	-----------------------------------------------------------------------------------
	'''
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
		structure = Bio.PDB.PDBParser().get_structure('temp1', 'temp1.pdb')
		dssp = Bio.PDB.DSSP(structure[0], 'temp1.pdb', acc_array='Wilke')
		ppb = Bio.PDB.Polypeptide.PPBuilder()
		chain = ppb.build_peptides(structure, aa_only=False)[0]
		SS = []
		for aa in dssp:
			if aa[2] == 'G' or aa[2] == 'H' or aa[2] == 'I': SSname = 'H'
			elif aa[2] == 'B' or aa[2] == 'E': SSname = 'S'
			else: SSname = 'L'
			SS.append(SSname)
		# Adjust End
		try:
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
			pose.dump_pdb('backbone.pdb')
			os.remove('temp2.pdb')
		except:
			os.remove('temp1.pdb')
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
			try:
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
				if Filter('backbone.pdb'):
					break
				else:
					os.remove('backbone.pdb')
			except:
				continue

def main():
	if args.train:
		LSTM_GAN('train')
	elif args.fragments:
		LSTM_GAN('predict')
		RD = RosettaDesign()
		RD.flxbb('backbone.pdb')
		Fragments('structure.pdb', sys.argv[1])		
	else:
		LSTM_GAN('predict')
		RD = RosettaDesign()
		RD.flxbb('backbone.pdb')

if __name__ == '__main__': main()
