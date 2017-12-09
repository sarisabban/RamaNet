import os , re , time , datetime , random , requests , urllib.request , bs4 , Bio.PDB
from pyrosetta import *
from pyrosetta.toolbox import *
init()

def Draw(SecondaryStructureString):
	''' Draws a protein topology given its secondary structure and distances '''
	''' Generates the DeNovo.pdb file '''
	#Length of structure
	length = len(SecondaryStructureString)
	#Construct primary structure made of Valines
	Val = str()
	for itr in range(length):
		itr = 'V'
		Val = Val + itr
	pose = pose_from_sequence(Val)
	for repeat in range(1):
		#Apply torsion angles
		count = 0
		for resi in SecondaryStructureString:
			count += 1
			if resi == 'H':
				pose.set_phi(int(count) , -57.8)	#From http://www.cryst.bbk.ac.uk/PPS2/course/section8/ss-960531_5.html
				pose.set_psi(int(count) , -47.0)	#From http://www.cryst.bbk.ac.uk/PPS2/course/section8/ss-960531_5.html
			elif resi == 'S':
				pose.set_phi(int(count) , -120)		#From http://www.cryst.bbk.ac.uk/PPS2/course/section8/ss-960531_10.html#HEADING9
				pose.set_psi(int(count) , 120)		#From http://www.cryst.bbk.ac.uk/PPS2/course/section8/ss-960531_10.html#HEADING9
		#Add constraints option to pose
		constraints = pyrosetta.rosetta.protocols.simple_moves.ConstraintSetMover()
		constraints.constraint_file('constraints.cst')
		constraints.add_constraints(True)
		constraints.apply(pose)
		#Score function with weight on only atom_pair_constraint
		scorefxn = ScoreFunction()
		scorefxn.set_weight(pyrosetta.rosetta.core.scoring.ScoreType.atom_pair_constraint , 	1.0)
#		scorefxn.set_weight(pyrosetta.rosetta.core.scoring.ScoreType.atom_pair_constraint , 	1.0)
#		scorefxn.set_weight(pyrosetta.rosetta.core.scoring.ScoreType.atom_pair_constraint , 	1.0)
#		scorefxn.set_weight(pyrosetta.rosetta.core.scoring.ScoreType.atom_pair_constraint , 	1.0)
		#Constraint relax to bring atoms together
		relax = pyrosetta.rosetta.protocols.relax.FastRelax()
		relax.set_scorefxn(scorefxn)
#		relax.constrain_relax_to_start_coords(True)
#		relax.constrain_coords(True)
		relax.apply(pose)
		#Normal FastRelax with constraints
#		scorefxn = pyrosetta.rosetta.core.scoring.ScoreFunctionFactory.create_score_function('ref2015_cst')
#		relax = pyrosetta.rosetta.protocols.relax.FastRelax()
#		relax.set_scorefxn(scorefxn)
#		relax.constrain_relax_to_start_coords(True)
#		relax.constrain_coords(True)
#		relax.show(pose)
#		relax.apply(pose)
		#Normal FastRelax without constraints
#		scorefxn = get_fa_scorefxn()
#		relax = pyrosetta.rosetta.protocols.relax.FastRelax()
#		relax.set_scorefxn(scorefxn)
#		relax.apply(pose)
	pose.dump_pdb('DeNovo.pdb')


#SS = 'LSSSSLLLLLLHHHHHHHHHHHHHHL' #AtomPair CA 1 CA 26 GAUSSIANFUNC 13.2 1.0	AtomPair CA 6 CA 11 GAUSSIANFUNC 9.7 1.0
#SS = 'SSSSSSLHHHHHHHHHHHHHHHHL'
#SS = 'SSSSSSLHHHL'
#SS = 'SSSSSSLLL'
#Draw(SS)
#distance i. 1 and n. CA, i. 26 and n. CA
#distance i. 6 and n. CA, i. 11 and n. CA

SS = 'LSSSSLLLLLLHHHHHHHHHHHHHHLSSSSSSLHHHHHHHHHHHHHHHHLSSSSSSLHHHLSSSSSSLLL'
Draw(SS)
