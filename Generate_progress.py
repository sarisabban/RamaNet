#!/usr/bin/python3

import os
import sys
from pyrosetta import *
from pyrosetta.toolbox import *
init()

def findmaxCST(filename):
	data = open(filename, 'r')
	cstALL = []
	for line in data:
		line = line.strip().split(';')
		cstALL.append(line[1::3])
	cst = []
	for value in cstALL:
		for item in value:
			try:
				item = float(item)
				cst.append(item)
			except:
				pass
	print(max(cst))
#findmaxCST('simple_dataset.csv')

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
	#Move amino acids angles
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
	#Score function with weight on only atom_pair_constraint
	scorefxnCST = ScoreFunction()
	scorefxnCST.set_weight(pyrosetta.rosetta.core.scoring.ScoreType.atom_pair_constraint, 1.0)
	#Constraint relax to bring atoms together
	relaxCST = pyrosetta.rosetta.protocols.relax.FastRelax()
	relaxCST.set_scorefxn(scorefxnCST)
	relaxCST.constrain_relax_to_start_coords(True)
	relaxCST.constrain_coords(True)
	#Normal FastRelax with constraints
	scorefxn = get_fa_scorefxn()
	relaxC = pyrosetta.rosetta.protocols.relax.FastRelax()
	relaxC.set_scorefxn(scorefxn)
	relaxC.constrain_relax_to_start_coords(True)
	relaxC.constrain_coords(True)
	#Normal FastRelax without constraints
	relax = pyrosetta.rosetta.protocols.relax.FastRelax()
	relax.set_scorefxn(scorefxn)
	#Run Relaxations
#	relaxCST.apply(pose)	#only brings structure together
	relaxC.apply(pose)	#Best on its own
#	relax.apply(pose)	#Best on its own
	pose.dump_pdb('Backbone.pdb')
	os.remove('constraints.cst')

def direct():
	for TheFile in os.listdir('results'):
		thefile = open('results/{}'.format(TheFile), 'r')
		#Remove zero ends
		newfile = open('results/X{}'.format(TheFile), 'a')
		for line in thefile:
			line = line.strip().split(';')
			newline = [round(float(line[0]), 3), round(float(line[1]), 3), round(float(line[2]), 3)]
			if newline[0] == 0.0 and newline[1] == 0.0 and newline[2] == 0.0:
				continue
			else:
				line = ';'.join(line)
				newfile.write(line + '\n')
		thefile.close()
		newfile.close()
		newfile = open('results/X{}'.format(TheFile), 'r')
		phiout = []
		psiout = []
		cstout = []
		for line in newfile:
			line = line.strip().split(';')
			phiout.append(float(line[0]))
			psiout.append(float(line[1]))
			cstout.append(float(line[2]))
		phiout = [x*360.0 for x in phiout]
		psiout = [x*360.0 for x in psiout]
		cstout = [x*54.505 for x in cstout]#small database 54.505 (old 207.801)
		data = (phiout, psiout, cstout)
		FoldPDB_PSC(data)
		Name = TheFile.split('.')[0]
		os.rename('Backbone.pdb', '{}.pdb'.format(TheFile))
	os.system('mv *.pdb results')
	os.system('rm results/*.txt')
	os.system('tar -zvcf results.tar.gz results')

def once():
	thefile = open('ori.txt', 'r')
	#Remove zero ends
	newfile = open('xx.txt', 'a')
	for line in thefile:
		line = line.strip().split(';')
		newline = [round(float(line[0]), 3), round(float(line[1]), 3), round(float(line[2]), 3)]
		if newline[0] == 0.0 and newline[1] == 0.0 and newline[2] == 0.0:
			continue
		else:
			line = ';'.join(line)
			newfile.write(line + '\n')
	thefile.close()
	newfile.close()
	newfile = open('xx.txt', 'r')
	phiout = []
	psiout = []
	cstout = []
	for line in newfile:
		line = line.strip().split(';')
		phiout.append(float(line[0]))
		psiout.append(float(line[1]))
		cstout.append(float(line[2]))
	phiout = [x*360.0 for x in phiout]
	psiout = [x*360.0 for x in psiout]
	cstout = [x*54.505 for x in cstout]#small database 54.505 (old 207.801)
	data = (phiout, psiout, cstout)
	FoldPDB_PSC(data)
	os.remove('xx.txt')

#direct()
once()
