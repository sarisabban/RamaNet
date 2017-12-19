import os , math , Bio.PDB ; from pyrosetta import * ; from pyrosetta.toolbox import * ; init()

def Draw(BPfile , CSTfile , RgCutoff):
	#Generate a starting structure
	temp = open('temp.pdb' , 'w')
	temp.write('ATOM      1  N   VAL A  1       25.945   4.358  33.648  1.00 22.51           N  \nATOM      2  CA  VAL A  1       26.375   4.305  35.016  1.00 24.82           C  \nATOM      3  C   VAL A  1       27.860   4.146  35.064  1.00 17.93           C  \nATOM      4  O   VAL A  1       28.451   3.503  34.206  1.00 27.99           O  \nATOM      5  CB  VAL A  1       25.647   3.121  35.836  1.00 38.86           C  \nATOM      6  CG1 VAL A  1       24.936   2.161  34.876  1.00 40.73           C  \nATOM      7  CG2 VAL A  1       26.659   2.335  36.692  1.00 39.90           C  ')
	temp.close()
	pose = pose_from_pdb('temp.pdb')
	for iteration in range(100):
		#Run the BluePrintBDR mover
		scorefxn = pyrosetta.rosetta.core.scoring.ScoreFunctionFactory.create_score_function('fldsgn_cen')
		scorefxn.set_weight(rosetta.core.scoring.atom_pair_constraint , 1.0)
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
		#Relax Structure
		scorefxn = get_fa_scorefxn()
		relax = pyrosetta.rosetta.protocols.relax.FastRelax()
		relax.set_scorefxn(scorefxn)
		relax.apply(pose)
		pose.dump_pdb('DeNovo.pdb')
		#Evaluation - Rg
		mass = list()
		Structure = open('DeNovo.pdb' , 'r')
		for line in Structure:
			line = line.split()
			try:
				if line[0] != 'ATOM':
					continue
				else:
					if line[-1] == 'C':
						mass.append(12.0107)
					elif line[-1] == 'O':
						mass.append(15.9994)
					elif line[-1] == 'N':
						mass.append(14.0067)
					elif line[-1] == 'H':
						mass.append(1.00794)
					else:
						continue
			except:
				continue
		coord = list()
		p = Bio.PDB.PDBParser()
		structure = p.get_structure('X', 'DeNovo.pdb')
		for model in structure:
			for chain in model:
				for residue in chain:
					for atom in residue:
						coord.append(atom.get_coord())
		xm = [(m * i , m * j , m * k) for (i , j , k) , m in zip(coord , mass)]
		tmass = sum(mass)
		rr = sum(mi * i + mj * j + mk * k for (i , j , k) , (mi , mj , mk) in zip(coord , xm))
		mm = sum((sum(i) / tmass) ** 2 for i in zip( * xm))
		rg = math.sqrt(rr / tmass - mm)
		print(rg)
		if rg <= RgCutoff:
			#Evaluation - Loops
			parser = Bio.PDB.PDBParser()
			structure = parser.get_structure('X' , 'DeNovo.pdb')
			model = structure[0]
			dssp = Bio.PDB.DSSP(model , 'DeNovo.pdb' , acc_array = 'Wilke')
			SS = list()
			for res in dssp:
				ss = res[2]
				if ss == '-' or ss == 'T' or ss == 'S':
					SS.append('L')
				elif ss == 'G' or ss == 'H' or ss == 'I' or ss == 'B' or ss == 'E':
					SS.append('NL')
			loop = SS.count('L')
			notloop = SS.count('NL')
			if loop >= notloop:
				os.remove('DeNovo.pdb')
				continue
			else:
				break
		else:
			os.remove('DeNovo.pdb')
			continue
	os.remove('temp.pdb')

#Draw('blueprint.bpf' , 'constraints.cst' , 15)
