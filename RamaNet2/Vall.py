import os

def Vall_to_CSV():
	''' Converts the vall.jul19.2011 database from Rosetta to a .csv file '''
	assert os.path.isfile('./vall.jul19.2011'),\
	'Make sure the vall.jul19.2011 file is in the same directory as this script'
	with open('vall.jul19.2011', 'r') as f:
		with open('vall.jul19.2011.csv', 'w') as F:
			F.write('PDBID,Chain,AA,SS,Sequence_Position,B_Factor,CA_x,CA_y,CA_z,CB_x,CB_y,CB_z,CEN_x,CEN_y,CEN_z,Phi,Psi,Omega,DSSP_phi,DSSP_psi,DSSP_sa,n_ali,pA,pC,pD,pE,pF,pG,pH,pI,pK,pL,pM,pN,pP,pQ,pR,pS,pT,pV,pW,pY,psA,psC,psD,psE,psF,psG,psH,psI,psK,psL,psM,psN,psP,psQ,psR,psS,psT,psV,psW,psY\n')
			for i in range(30): next(f)
			for line in f:
				line = line.strip().split()
				ID = line[0][:4]
				CH = line[0][-1]
				RE = ','.join(line[1:])
				line = '{},{},{}\n'.format(ID, CH, RE)
				F.write(line)

def Vall():
	'''
	Compile the PDB IDs, chains, phi, psi, omega, and SASA of all the structures
	from the Rosetta vall.jul19.2011 database into a .csv file
	 '''
	assert os.path.isfile('./vall.jul19.2011'),\
	'Make sure the vall.jul19.2011 file is in the same directory as this script'
	nx = 1490
	m  = 16800
	with open('vall.jul19.2011', 'r') as f:
		with open('Fragments_no_fill.csv', 'w') as F:
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

def Count(filename):
	''' Get the largest structure size for the .csv file header (nx value) '''
	sizes = []
	with open('dataset.csv', 'r') as f:
		for line in f:
			line = line.strip().split(',')
			sizes.append(len(line[2:])/6)
	print(max(sizes), len(sizes))

def Fill(filename):
    ''' Fills missing .csv table spaces with 0s '''
    with open(filename) as f:
	    with open('Fragments.csv', 'a') as F:
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

Vall()
Fill('Fragments.csv')
