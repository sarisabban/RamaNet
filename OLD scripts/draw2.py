import os , re , time , datetime , random , requests , urllib.request , bs4 , Bio.PDB , math
from pyrosetta import *
from pyrosetta.toolbox import *
init()

#Get torsion angels
parser = Bio.PDB.PDBParser()
structure = parser.get_structure('X' , 'test.pdb')
model = structure[0]
dssp = Bio.PDB.DSSP(model , 'test.pdb' , acc_array = 'Wilke')
phi = list()
psi = list()
for aa in dssp:
	phi.append(aa[4])
	psi.append(aa[5])
phi = [0.0 if x == 360.0 else x for x in phi]
psi = [0.0 if x == 360.0 else x for x in psi]
#Get secondary structures
parser = Bio.PDB.PDBParser()
structure = parser.get_structure('X' , 'test.pdb')
model = structure[0]
dssp = Bio.PDB.DSSP(model , 'test.pdb' , acc_array='Wilke')
SS = list()
for res in dssp:
	ss = res[2]
	if ss == '-' or ss == 'T' or ss == 'S':		#Loop (DSSP code is - or T or S)
		SS.append('L')
	elif ss == 'G' or ss == 'H' or ss == 'I':	#Helix (DSSP code is G or H or I)
		SS.append('H')
	elif ss == 'B' or ss == 'E':			#Sheet (DSSP code is B or E)
		SS.append('S')
SS = ''.join(SS)
'''
#Get Constraints
#https://www.rosettacommons.org/docs/latest/rosetta_basics/file_types/constraint-file
dssp = Bio.PDB.DSSP(model , 'test.pdb' , acc_array = 'Wilke')
for aa in dssp:
	length = aa[0]
structure = Bio.PDB.PDBParser(QUIET=True).get_structure('X' , 'test.pdb')
ppb = Bio.PDB.Polypeptide.PPBuilder()
Type = ppb.build_peptides(structure , aa_only=True)
model = Type
chain = model[0]
count = 1
for aa in range(length - 1):
	residue1 = chain[0]
	residue2 = chain[count]
	atom1 = residue1['CA']
	atom2 = residue2['CA']
	constraint = atom1 - atom2
	count += 1
	thefile = open('constraints.cst' , 'a')
	line = 'AtomPair CA 1 CA ' + str(count) +' GAUSSIANFUNC '+ str(constraint) +' 1.0\n'
	thefile.write(line)
	thefile.close()
'''
def Draw1(SecondaryStructureString , phi_list , psi_list):
	''' Draws a protein topology given its secondary structure '''
	''' Generates the DeNovo.pdb file '''
	#Length of structure
	length = len(SecondaryStructureString)
	#Construct primary structure made of Valines
	Val = str()
	for itr in range(length):
		itr = 'V'
		Val = Val + itr
	pose = pose_from_sequence(Val)
	#Apply torsion angles
	count = 0
	for resi , p , s in zip(SecondaryStructureString , phi_list , psi_list):
		count += 1
		pose.set_phi(int(count) , p)
		pose.set_psi(int(count) , s)
	#Add constraints option to pose
	constraints = pyrosetta.rosetta.protocols.simple_moves.ConstraintSetMover()
	constraints.constraint_file('constraints.cst')
	constraints.add_constraints(True)
	constraints.apply(pose)
	#Score function with weight on only atom_pair_constraint
	scorefxn = ScoreFunction()
	scorefxn.set_weight(pyrosetta.rosetta.core.scoring.ScoreType.atom_pair_constraint , 	1.0)
#	scorefxn.set_weight(pyrosetta.rosetta.core.scoring.ScoreType.atom_pair_constraint , 	1.0)
#	scorefxn.set_weight(pyrosetta.rosetta.core.scoring.ScoreType.atom_pair_constraint , 	1.0)
#	scorefxn.set_weight(pyrosetta.rosetta.core.scoring.ScoreType.atom_pair_constraint , 	1.0)
	#Constraint relax to bring atoms together
	relax = pyrosetta.rosetta.protocols.relax.FastRelax()
	relax.set_scorefxn(scorefxn)
#	relax.constrain_relax_to_start_coords(True)
#	relax.constrain_coords(True)
	relax.apply(pose)
	#Normal FastRelax with constraints
#	scorefxn = pyrosetta.rosetta.core.scoring.ScoreFunctionFactory.create_score_function('ref2015_cst')
#	relax = pyrosetta.rosetta.protocols.relax.FastRelax()
#	relax.set_scorefxn(scorefxn)
#	relax.constrain_relax_to_start_coords(True)
#	relax.constrain_coords(True)
#	relax.show(pose)
#	relax.apply(pose)
	#Normal FastRelax without constraints
#	scorefxn = get_fa_scorefxn()
#	relax = pyrosetta.rosetta.protocols.relax.FastRelax()
#	relax.set_scorefxn(scorefxn)
#	relax.apply(pose)
	pose.dump_pdb('DeNovo.pdb')

def Draw2(SS , phi , psi):
#	https://www.rosettacommons.org/docs/latest/scripting_documentation/RosettaScripts/Movers/movers_pages/BluePrintBDRMover
	#5 - Run the BluePrintBDR mover
	mover = pyrosetta.rosetta.protocols.fldsgn.BluePrintBDR()
	mover.num_fragpick(200)
	mover.use_fullmer(True)
	mover.use_abego_bias(True)
	mover.use_sequence_bias(False)
	mover.max_linear_chainbreak(0.07)
	mover.ss_from_blueprint(True)
	mover.dump_pdb_when_fail('')
	mover.set_constraints_NtoC(-1.0)
	mover.set_blueprint('blueprint')
	mover.apply(pose)
	os.remove('blueprint')
#	mover.set_constraint_file('structure.cst')
#	mover.scorefunction('ref2015_cst')
	PoseFinal.dump_pdb('DeNovo.pdb')
	os.remove('temp.pdb')




Draw1(SS , phi , psi)







