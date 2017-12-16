import os ; from pyrosetta import * ; from pyrosetta.toolbox import * ; init()

def Draw(BPfile , CSTfile):
	#Generate a starting structure
	temp = open('temp.pdb' , 'w')
	temp.write('ATOM      1  N   VAL A  1       25.945   4.358  33.648  1.00 22.51           N  \nATOM      2  CA  VAL A  1       26.375   4.305  35.016  1.00 24.82           C  \nATOM      3  C   VAL A  1       27.860   4.146  35.064  1.00 17.93           C  \nATOM      4  O   VAL A  1       28.451   3.503  34.206  1.00 27.99           O  \nATOM      5  CB  VAL A  1       25.647   3.121  35.836  1.00 38.86           C  \nATOM      6  CG1 VAL A  1       24.936   2.161  34.876  1.00 40.73           C  \nATOM      7  CG2 VAL A  1       26.659   2.335  36.692  1.00 39.90           C  ')
	temp.close()
	pose = pose_from_pdb('temp.pdb')
	os.remove('temp.pdb')
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


Draw('blueprint.bpf' , 'constraints.cst')
