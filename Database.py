#!/usr/bin/python3

import os , math , gzip , Bio.PDB , Bio.pairwise2 , tqdm

def Database(TempDIR , FinalDIR):
	''' Downloads the entire PDB database from https://www.wwpdb.org/, moves all files into one directory, then uncompresses all the files '''
	''' Generates a directory which contains all .PDB structure files '''
	os.system('rsync -rlpt -v -z --delete --port=33444 rsync.wwpdb.org::ftp/data/structures/divided/pdb/ ./' + TempDIR)
	os.mkdir(FinalDIR)
	filelist = os.listdir(TempDIR)
	print('\x1b[32m' + 'Download complete' + '\x1b[0m')
	print('\x1b[32m' + 'Moving files' + '\x1b[0m')
	for directories in tqdm.tqdm(filelist):
		files = os.listdir(TempDIR + '/' + directories)
		for afile in files:
			location = (TempDIR + '/' + directories + '/' + afile)
			os.rename(location , FinalDIR + '/' + afile)
	os.system('rm -r ./' + TempDIR)

def Extract(directory):
	''' Extracts all the .ent.gz files and separate all chains and save them into seperate .pdb files '''
	''' Replaces each .ent.gz file with the .pdb file of each chain '''
	current = os.getcwd()
	pdbfilelist = os.listdir(directory)
	io = Bio.PDB.PDBIO()
	os.chdir(directory)
	print('\x1b[32m' + 'Extracting files' + '\x1b[0m')
	for TheFile in tqdm.tqdm(pdbfilelist):
		try:
			TheName = TheFile.split('.')[0].split('pdb')[1].upper()						#Open file
			InFile = gzip.open(TheFile, 'rt')											#Extract file
			structure = Bio.PDB.PDBParser(QUIET = True).get_structure(TheName , InFile)	#Separate chains and save to different files
			count = 0
			for chain in structure.get_chains():
				io.set_structure(chain)
				io.save(structure.get_id() + '_' + chain.get_id() + '.pdb')
			os.remove(TheFile)
		except Exception as TheError:
			print('\x1b[31m' + '[-] Failed to extract' + '\t' + thefile.upper() , '\x1b[33m' + str(TheError) + '\x1b[0m')
			os.remove(TheFile)
	os.chdir(current)

def NonProtein(directory):
	''' Remove non-protein structures '''
	current = os.getcwd()
	pdbfilelist = os.listdir(directory)
	os.chdir(directory)
	print('\x1b[32m' + 'Deleting none-protein structures' + '\x1b[0m')
	for TheFile in tqdm.tqdm(pdbfilelist):
		structure = Bio.PDB.PDBParser(QUIET = True).get_structure('X' , TheFile)
		ppb = Bio.PDB.Polypeptide.PPBuilder()
		Type = ppb.build_peptides(structure , aa_only = True)
		if Type == []:																	#Non-protein structures have Type = []
			os.remove(TheFile)
		else:
			continue
	os.chdir(current)

def Size(directory , Size_From , Size_To):
	''' Remove 80AA < structures < 150AA '''
	current = os.getcwd()
	pdbfilelist = os.listdir(directory)
	os.chdir(directory)
	print('\x1b[32m' + 'Removing structure sizes less than 80 amino acids or larger than 150 amino acids' + '\x1b[0m')
	for TheFile in tqdm.tqdm(pdbfilelist):
		try:
			parser = Bio.PDB.PDBParser()
			structure = parser.get_structure('X' , TheFile)
			model = structure[0]
			dssp = Bio.PDB.DSSP(model , TheFile , acc_array = 'Wilke')
			for aa in dssp:																#Identify final structure's length
				length = aa[0]
			if length >= int(Size_To) or length <= int(Size_From):
				os.remove(TheFile)
		except:
			print('\x1b[31m' + 'Error in finding protein size' + '\x1b[0m')
	os.chdir(current)

def Break(directory):
	''' Remove structures with a broken (non-continuous) chains '''
	current = os.getcwd()
	pdbfilelist = os.listdir(directory)
	os.chdir(directory)
	print('\x1b[32m' + 'Removing structures with non-continuous chains' + '\x1b[0m')
	for TheFile in tqdm.tqdm(pdbfilelist):
		structure = Bio.PDB.PDBParser(QUIET = True).get_structure('X' , TheFile)
		ppb = Bio.PDB.Polypeptide.PPBuilder()
		Type = ppb.build_peptides(structure , aa_only = True)
		try:
			x = Type[1]
			os.remove(TheFile)
		except:
			continue
	os.chdir(current)

def Loops(directory , LoopLength):
	''' Remove structures that have loops that are larger than a spesific length '''
	current = os.getcwd()
	pdbfilelist = os.listdir(directory)
	os.chdir(directory)
	print('\x1b[32m' + 'Removing structures with long loops' + '\x1b[0m')
	for TheFile in tqdm.tqdm(pdbfilelist):
		try:
			parser = Bio.PDB.PDBParser()
			structure = parser.get_structure('X' , TheFile)
			model = structure[0]
			dssp = Bio.PDB.DSSP(model , TheFile , acc_array = 'Wilke')
			SS = list()
			for res in dssp:
				ss = res[2]
				if ss == '-' or ss == 'T' or ss == 'S':										#Loop (DSSP code is - or T or S)
					SS.append('L')
				else:
					SS.append('.')
			loops = ''.join(SS).split('.')
			loops = [item for item in loops if item] 
			LargeLoop = None
			for item in loops:
				if len(item) <= LoopLength:
					continue
				else:
					LargeLoop = 'LargeLoop'
			if LargeLoop == 'LargeLoop':
				os.remove(TheFile)
			else:
				continue
		except:
			os.remove(TheFile)
	os.chdir(current)

def Renumber(directory):
	''' Renumber structures starting at 1 '''
	current = os.getcwd()
	pdbfilelist = os.listdir(directory)
	os.chdir(directory)
	print('\x1b[32m' + 'Renumbering structures' + '\x1b[0m')
	for TheFile in tqdm.tqdm(pdbfilelist):
		pdb = open(TheFile , 'r')
		PDB = open(TheFile + 'X' , 'w')
		count = 0
		num = 0
		AA2 = None
		for line in pdb:
			count += 1																	#Sequencially number atoms
			AA1 = line[23:27]															#Sequencially number residues
			if not AA1 == AA2:
				num += 1			
			final_line = line[:7] + '{:4d}'.format(count) + line[11:17] + line[17:21] + 'A' + '{:4d}'.format(num) + line[26:]	#Update each line to have its atoms and residues sequencially labeled, as well as being in chain A
			AA2 = AA1
			PDB.write(final_line)														#Write to new file called motif.pdb
		PDB.close()
		os.remove(TheFile)
		os.rename(TheFile + 'X' , TheFile)
	os.chdir(current)

def RMSD(directory , RMSDcutoff):
	''' Remove structures that are similar to each other '''
	current = os.getcwd()
	pdbfilelist = os.listdir(directory)
	os.chdir(directory)
	print('\x1b[32m' + 'Removing structure with similar RMSD' + '\x1b[0m')
	for File1 in tqdm.tqdm(pdbfilelist):
		for File2 in pdbfilelist:
			if File1 == File2:
				continue
			else:
				try:
					#First structure
					type1 = Bio.PDB.Polypeptide.PPBuilder().build_peptides(Bio.PDB.PDBParser(QUIET = True).get_structure('X' , File1) , aa_only = True)
					length1 = type1[-1][-1].get_full_id()[3][1]
					fixed = [atom['CA'] for atom in type1[0]]
					#Second structure
					type2 = Bio.PDB.Polypeptide.PPBuilder().build_peptides(Bio.PDB.PDBParser(QUIET = True).get_structure('X' , File2) , aa_only = True)
					length2 = type2[-1][-1].get_full_id()[3][1]
					moving = [atom['CA'] for atom in type2[0]]
					#Choose the length of the smallest structure
					lengths = [length1 , length2]
					smallest = min(int(item) for item in lengths)
					#Find RMSD
					sup = Bio.PDB.Superimposer()
					sup.set_atoms(fixed[:smallest] , moving[:smallest])
					sup.apply(Bio.PDB.PDBParser(QUIET = True).get_structure('X' , File2)[0].get_atoms())
					RMSD = round(sup.rms , 4)
					print(File1 , File2 , RMSD)
					#Delete similar structures
					if RMSD < RMSDcutoff:
						os.remove(File2)
				except:
					continue
	os.chdir(current)

def Sequence(directory , Cutoff):
	''' Remove structures that have similar sequences, which means they most likely have similar structures '''
	current = os.getcwd()
	pdbfilelist = os.listdir(directory)
	os.chdir(directory)
	print('\x1b[32m' + 'Measuring sequence similarity' + '\x1b[0m')
	for File1 in tqdm.tqdm(pdbfilelist):
			for File2 in pdbfilelist:
				try:
					if File1 == File2:
						continue
					else:
						seq1 = Bio.PDB.Polypeptide.PPBuilder().build_peptides(Bio.PDB.PDBParser(QUIET = True).get_structure('X' , File1) , aa_only = True)[0].get_sequence()
						seq2 = Bio.PDB.Polypeptide.PPBuilder().build_peptides(Bio.PDB.PDBParser(QUIET = True).get_structure('X' , File2) , aa_only = True)[0].get_sequence()
						alignment = Bio.pairwise2.align.globalxx(seq1 , seq2)
						total = alignment[0][4]
						similarity = alignment[0][2]
						percentage = (similarity * 100) / total
						if percentage > Cutoff:
							os.remove(File2)
				except:
					continue
	os.chdir(current)

def Rg(directory , RGcutoff):
	''' Remove structures that are below the Raduis of Gyration's value '''
	current = os.getcwd()
	pdbfilelist = os.listdir(directory)
	os.chdir(directory)
	print('\x1b[32m' + 'Removing structure low Rg values' + '\x1b[0m')
	for TheFile in tqdm.tqdm(pdbfilelist):
		mass = list()
		Structure = open(TheFile , 'r')
		for line in Structure:
			line = line.split()
			if line[0] == 'TER' or line[0] == 'END':
				continue
			else:
				if line[-1] == 'C':
					mass.append(12.0107)
				elif line[-1] == 'O':
					mass.append(15.9994)
				elif line[-1] == 'N':
					mass.append(14.0067)
				elif line[-1] == 'S':
					mass.append(32.0650)
				elif line[-1] == 'H':
					mass.append(1.00794)
				else:
					continue
		coord = list()
		p = Bio.PDB.PDBParser()
		structure = p.get_structure('X', TheFile)
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
		if rg <= RGcutoff:
			os.remove(TheFile)
		else:
			continue
	os.chdir(current)

def SS(directory):
	''' Get the secondary structures '''
	''' Generates a list with each amino acid's secondary strucure for each protein in a directory '''
	current = os.getcwd()
	pdbfilelist = os.listdir(directory)
	os.chdir(directory)
	print('\x1b[32m' + 'Getting the secondary structures of each protein' + '\x1b[0m')
	count = 1
	for TheFile in tqdm.tqdm(pdbfilelist):
		parser = Bio.PDB.PDBParser()
		structure = parser.get_structure('X' , TheFile)
		model = structure[0]
		dssp = Bio.PDB.DSSP(model , TheFile , acc_array='Wilke')
		SS = list()
		for res in dssp:
			ss = res[2]
			if ss == '-' or ss == 'T' or ss == 'S':										#Loop (DSSP code is - or T or S)
				SS.append('L')
			elif ss == 'G' or ss == 'H' or ss == 'I':									#Helix (DSSP code is G or H or I)
				SS.append('H')
			elif ss == 'B' or ss == 'E':												#Sheet (DSSP code is B or E)
				SS.append('S')
		SS = ['1' if x == 'L' else x for x in SS]
		SS = ['2' if x == 'H' else x for x in SS]
		SS = ['3' if x == 'S' else x for x in SS]
		addition = 150 - len(SS)
		for zeros in range(addition):
			SS.append('0')
		line = ';'.join(SS)
		data = open('SS' , 'a')
		data.write(str(count) + ';' + TheFile.split('.')[0] + ';' + line + '\n')
		data.close()
		count += 1
	os.chdir(current)
	os.rename(directory + '/SS' , 'SS')

def Distances(directory):
	''' Measures distances between the first amino acid and all the others '''
	''' Generates a list with the distances '''
	current = os.getcwd()
	pdbfilelist = os.listdir(directory)
	os.chdir(directory)
	print('\x1b[32m' + 'Measuring distances' + '\x1b[0m')
	for TheFile in tqdm.tqdm(pdbfilelist):
		parser = Bio.PDB.PDBParser()
		structure = parser.get_structure('X' , TheFile)
		model = structure[0]
		dssp = Bio.PDB.DSSP(model , TheFile , acc_array = 'Wilke')
		for aa in dssp:																	#Identify final structure's length
			length = aa[0]
		structure = Bio.PDB.PDBParser(QUIET = True).get_structure('X' , TheFile)
		ppb = Bio.PDB.Polypeptide.PPBuilder()
		Type = ppb.build_peptides(structure , aa_only = True)
		model = Type
		chain = model[0]
		positions = list(range(length // 10, length, length // 10))  
		distances = list()
		for res in positions:
			try:
				residue1 = chain[0]
				residue2 = chain[res + 1]
				atom1 = residue1['CA']
				atom2 = residue2['CA']
				distance = atom1-atom2
				distances.append(str(distance))
			except:
				continue

		line = ';'.join(distances)
		data = open('DI' , 'a')
		data.write(line + '\n')
		data.close()
	os.chdir(current)
	os.rename(directory + '/DI' , 'DI')

def PutTogether(SS , DI):
	sec = open(SS , 'r')
	dis = open(DI , 'r')


#---------------------------------------------------------------------------------------------------------------------------------------
#Protocol to isolate specific types of structures
#Database('DATABASE' , 'PDBDatabase')	# 1. Download the PDB database
#Extract('PDBDatabase')			# 2. Extract files
#NonProtein('PDBDatabase')		# 3. Remove non-protein structures
#Size('PDBDatabase' , 80 , 150)		# 4. Remove structures less than or larger than a specified amino acid leangth
#Break('PDBDatabase')			# 5. Remove structure with broken chains
#Loops('PDBDatabase' , 10)		# 6. Remove structures that have loops that are larger than a spesific length
#Renumber('PDBDatabase')		# 7. Renumber structures starting at amino acid 1
#Sequence('PDBDatabase' , 75)		# 8. Align the sequences of each structure to each structure, remove structures with similar sequences that fall above a user defined percentage
#RMSD('PDBDatabase' , 5)		# 9. Measure RMSD of each structure to each structure, remove if RMSD < specified value (CODE IS NOT VERY RELIABLE)
#Rg('PDBDatabase' , 15)			# 10. Remove structures that are below a specified Raduis of Gyration value

#Protocol to extract specific information from isolated structures
#SS('PDBDatabase')			# 1. Get the secondary structures
#Distances('PDBDatabase')		# 2. Measure distances between the first amino acid and all the others
PutTogether('SS' , 'DI')
