#!/usr/bin/python

import os
import sys
import Bio.PDB
from pyrosetta import *
from pyrosetta.toolbox import *
init()

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
	pose.dump_pdb('temp.pdb')
	structure = Bio.PDB.PDBParser().get_structure('temp', 'temp.pdb')
	dssp = Bio.PDB.DSSP(structure[0], 'temp.pdb', acc_array='Wilke')
	ppb = Bio.PDB.Polypeptide.PPBuilder()
	chain = ppb.build_peptides(structure, aa_only=False)[0]
	SS = []
	for aa in dssp:
		if aa[2] == 'G' or aa[2] == 'H' or aa[2] == 'I':
			SSname = 'H'
		elif aa[2] == 'B' or aa[2] == 'E':
			SSname = 'S'
		else:
			SSname = 'L'
		SS.append(SSname)
	# Adjust End
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
	os.remove('temp.pdb')
	pose = pose_from_pdb('temp2.pdb')
	relax.apply(pose)
	pose.dump_pdb('Backbone.pdb')
	os.remove('temp2.pdb')

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
#	if Hnum+Snum < Lnum:
#		choice = False
#		print('Helix+Sheet:', Hnum+Snum, 'Loop:', Lnum)
	# Secondary structures number filter
	H = SS.count('H')
	S = SS.count('S')
	L = SS.count('L')
	if H+S < L:
		choice = False
#		print('Helix+Sheet Residues:', H+S, 'Loop residues', L)
	# SASA filter
	Surface = SASA.count('S')
	Boundary = SASA.count('B')
	Core = SASA.count('C')
	percent = (Core*100)/(Surface+Boundary+Core)
	if percent < 15:
		choice = False
#		print('SASA Core Percentage:', percent)
	# CST filter
	MaxCST = max(CST)
	if MaxCST > 88:
		choice = False
#		print('Max Constraint:', MaxCST)
	return(choice)

def main(directory):
	for TheFile in os.listdir(directory):
		newfile = open('{}/{}'.format(directory, TheFile), 'r')
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
		Name = TheFile.split('.')[0]
		os.rename('Backbone.pdb', '{}.pdb'.format(TheFile))
		if Filter('{}.pdb'.format(TheFile)):
			os.system('mv {}.pdb {}'.format(TheFile, 'OK'))
	os.system('mv *.pdb {}'.format(directory))
	os.mkdir('{}/good'.format(directory))
	os.mkdir('{}/bad'.format(directory))
	os.system('rm {}/*.txt'.format(directory))

if __name__ == '__main__': main(sys.argv[1])
