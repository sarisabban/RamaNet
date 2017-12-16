import os
from pyrosetta import *
from pyrosetta.toolbox import *
init()
'''
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
#Blueprint file
parser = Bio.PDB.PDBParser()
structure = parser.get_structure('X' , 'test.pdb')
model = structure[0]
dssp = Bio.PDB.DSSP(model , 'test.pdb' , acc_array = 'Wilke')
count = 1
for res in dssp:
	ss = res[2]
	if ss == '-' or ss == 'T' or ss == 'S':		#Loop (DSSP code is - or T or S)
		thefile = open('blueprint.bpf' , 'a')
		thefile.write('0 V LX R\n')
		thefile.close()
	elif ss == 'G' or ss == 'H' or ss == 'I':	#Helix (DSSP code is G or H or I)
		thefile = open('blueprint.bpf' , 'a')
		thefile.write('0 V HX R\n')
		thefile.close()
	elif ss == 'B' or ss == 'E':			#Sheet (DSSP code is B or E)
		thefile = open('blueprint.bpf' , 'a')
		thefile.write('0 V EX R\n')
		thefile.close()
	count += 1
with open('blueprint.bpf', 'r') as fin:
	data = fin.read().splitlines(True)
with open('blueprint.bpf', 'w') as fout:
	fout.writelines(data[1:])
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
	line = 'AtomPair CA 1 CA ' + str(count) +' GAUSSIANFUNC '+ str(constraint) +' 2.0 TAG\n'
	thefile.write(line)
	thefile.close()
'''
def Draw(BPfile , CSTfile):
	#https://www.rosettacommons.org/docs/latest/scripting_documentation/RosettaScripts/Movers/movers_pages/BluePrintBDRMover
	#Generate a starting structure
	temp = open('temp.pdb' , 'w')
	temp.write('ATOM      1  N   VAL A  1       25.945   4.358  33.648  1.00 22.51           N  \nATOM      2  CA  VAL A  1       26.375   4.305  35.016  1.00 24.82           C  \nATOM      3  C   VAL A  1       27.860   4.146  35.064  1.00 17.93           C  \nATOM      4  O   VAL A  1       28.451   3.503  34.206  1.00 27.99           O  \nATOM      5  CB  VAL A  1       25.647   3.121  35.836  1.00 38.86           C  \nATOM      6  CG1 VAL A  1       24.936   2.161  34.876  1.00 40.73           C  \nATOM      7  CG2 VAL A  1       26.659   2.335  36.692  1.00 39.90           C  ')
	temp.close()
	pose = pose_from_pdb('temp.pdb')
	#Run the BluePrintBDR mover
	scorefxn = pyrosetta.rosetta.core.scoring.ScoreFunctionFactory.create_score_function('ref2015_cst')
	mover = pyrosetta.rosetta.protocols.fldsgn.BluePrintBDR()
	mover.num_fragpick(200)
	mover.use_fullmer(True)
	mover.use_abego_bias(True)
	mover.use_sequence_bias(False)
	mover.max_linear_chainbreak(0.07)
	mover.ss_from_blueprint(True)
	mover.dump_pdb_when_fail('')
	mover.set_constraints_NtoC(-1.0)
	mover.set_blueprint(BPfile)
	mover.set_constraint_file(CSTfile)
	mover.scorefunction(scorefxn)
	mover.apply(pose)
	pose.dump_pdb('DeNovo.pdb')
	os.remove('temp.pdb')




Draw('blueprint.bpf' , 'constraints.cst')


















def Draw2(SecondaryStructureString , CSTfile):
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
	#Fold topology
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
	constraints.constraint_file(CSTfile)
	constraints.add_constraints(True)
	constraints.apply(pose)
	#Score function with weight on only atom_pair_constraint
	scorefxn = ScoreFunction()
	scorefxn.set_weight(pyrosetta.rosetta.core.scoring.ScoreType.atom_pair_constraint , 	1.0)
	#Constraint relax to bring atoms together
	relax = pyrosetta.rosetta.protocols.relax.FastRelax()
	relax.set_scorefxn(scorefxn)
	relax.constrain_relax_to_start_coords(True)
	relax.constrain_coords(True)
	relax.apply(pose)
	#Normal FastRelax with constraints
	scorefxn = pyrosetta.rosetta.core.scoring.ScoreFunctionFactory.create_score_function('ref2015_cst')
	relax = pyrosetta.rosetta.protocols.relax.FastRelax()
	relax.set_scorefxn(scorefxn)
	relax.constrain_relax_to_start_coords(True)
	relax.constrain_coords(True)
	relax.apply(pose)
	#Normal FastRelax without constraints
	scorefxn = get_fa_scorefxn()
	relax = pyrosetta.rosetta.protocols.relax.FastRelax()
	relax.set_scorefxn(scorefxn)
	relax.apply(pose)
	#Export
	pose.dump_pdb('DeNovo.pdb')
SS = 'LSSSSLLLLLLLHHHHHHHHHHHLLHHHLLSSSSSSLLHHHHHHHHHHHHHHLLSSSSSSLHHHLSSSSSSLLL'
#Draw2(SS , 'constraints.cst')

