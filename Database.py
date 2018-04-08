#!/usr/bin/python3

import os , math , gzip , Bio.PDB , Bio.pairwise2 , tqdm
#from pyrosetta import *
#from pyrosetta.toolbox import *
#init()

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
	print('\x1b[32m' + 'Measuring sequence similarity, round 1/2' + '\x1b[0m')
	for File1 in tqdm.tqdm(pdbfilelist):
		for File2 in pdbfilelist:
			try:
				if File1 == File2:
					continue
				else:
					if File1.split('.')[0].split('_')[0] == File2.split('.')[0].split('_')[0]:
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
	print('\x1b[32m' + 'Measuring sequence similarity, round 2/2' + '\x1b[0m')
	for File1 in tqdm.tqdm(pdbfilelist):
		for File2 in pdbfilelist:
			try:
				if File1 == File2:
					continue
				else:
					if File1.split('.')[0].split('_')[0][:3] == File2.split('.')[0].split('_')[0][:3]:
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

def DatasetR(directory):
	''' Get the secondary structures and distances '''
	''' Generates a the dataR.csv with each amino acid's secondary strucure and 10 distances between the first amino acid's CA atom and others for each protein in a directory '''
	current = os.getcwd()
	pdbfilelist = os.listdir(directory)
	os.chdir(directory)
	print('\x1b[32m' + 'Getting the secondary structures of each protein' + '\x1b[0m')
	data = open('dataR.csv' , 'a')
	data.write(';PDB_ID;1;2;3;4;5;6;7;8;9;10;11;12;13;14;15;16;17;18;19;20;21;22;23;24;25;26;27;28;29;30;31;32;33;34;35;36;37;38;39;40;41;42;43;44;45;46;47;48;49;50;51;52;53;54;55;56;57;58;59;60;61;62;63;64;65;66;67;68;69;70;71;72;73;74;75;76;77;78;79;80;81;82;83;84;85;86;87;88;89;90;91;92;93;94;95;96;97;98;99;100;101;102;103;104;105;106;107;108;109;110;111;112;113;114;115;116;117;118;119;120;121;122;123;124;125;126;127;128;129;130;131;132;133;134;135;136;137;138;139;140;141;142;143;144;145;146;147;148;149;150;Distance_1;Distance_2;Distance_3;Distance_4;Distance_5;Distance_6;Distance_7;Distance_8;Distance_9;Distance_10\n')
	data.close()
	count = 1
	for TheFile in tqdm.tqdm(pdbfilelist):
		try:
			structure = Bio.PDB.PDBParser().get_structure('X' , TheFile)
			model = structure[0]
			dssp = Bio.PDB.DSSP(model , TheFile , acc_array = 'Wilke')
			length = [aa[0] for aa in dssp][-1]			#Identify final structure's length
			SS = list()
			for res in dssp:
				ss = res[2]
				if ss == '-' or ss == 'T' or ss == 'S':		#Loop (DSSP code is - or T or S)
					SS.append('L')
				elif ss == 'G' or ss == 'H' or ss == 'I':	#Helix (DSSP code is G or H or I)
					SS.append('H')
				elif ss == 'B' or ss == 'E':			#Sheet (DSSP code is B or E)
					SS.append('S')
			SS = ['1' if x == 'L' else x for x in SS]
			SS = ['2' if x == 'H' else x for x in SS]
			SS = ['3' if x == 'S' else x for x in SS]
			addition = 150 - len(SS)
			for zeros in range(addition):
				SS.append('0')
			SSline =  ';'.join(SS)
			chain = Bio.PDB.Polypeptide.PPBuilder().build_peptides(structure , aa_only = True)[0]
			positions = [(i+1)*(length//10) for i in range(10)]
			distances = list()
			for res in positions:
				try:
					residue1 = chain[0]
					residue2 = chain[res - 1]
					atom1 = residue1['CA']
					atom2 = residue2['CA']
					distance = atom1-atom2
					distances.append(str(distance))
				except:
					continue
			if distances == []:
				continue
			elif len(distances) != 10:
				continue
			DIline = ';'.join(distances)
			data = open('dataR.csv' , 'a')
			data.write(str(count) + ';' + TheFile.split('.')[0] + ';' + SSline + ';' + DIline + '\n')
			data.close()
			count += 1
		except:
			continue
	os.chdir(current)
	os.rename(directory + '/dataR.csv' , 'dataR.csv')

def DatasetCA(directory):
	''' Get each residue's CA atom's XYZ coordinates '''
	''' Generates a the dataCA.csv with the XYZ coordinates of the CA atom for each amino acid '''
	current = os.getcwd()
	pdbfilelist = os.listdir(directory)
	os.chdir(directory)
	print('\x1b[32m' + "Getting the CA atom's XYZ coordinates" + '\x1b[0m')
	data = open('dataCA.csv' , 'a')
	data.write(';PDB_ID;X_1;Y_1;Z_1;X_2;Y_2;Z_2;X_3;Y_3;Z_3;X_4;Y_4;Z_4;X_5;Y_5;Z_5;X_6;Y_6;Z_6;X_7;Y_7;Z_7;X_8;Y_8;Z_8;X_9;Y_9;Z_9;X_10;Y_10;Z_10;X_11;Y_11;Z_11;X_12;Y_12;Z_12;X_13;Y_13;Z_13;X_14;Y_14;Z_14;X_15;Y_15;Z_15;X_16;Y_16;Z_16;X_17;Y_17;Z_17;X_18;Y_18;Z_18;X_19;Y_19;Z_19;X_20;Y_20;Z_20;X_21;Y_21;Z_21;X_22;Y_22;Z_22;X_23;Y_23;Z_23;X_24;Y_24;Z_24;X_25;Y_25;Z_25;X_26;Y_26;Z_26;X_27;Y_27;Z_27;X_28;Y_28;Z_28;X_29;Y_29;Z_29;X_30;Y_30;Z_30;X_31;Y_31;Z_31;X_32;Y_32;Z_32;X_33;Y_33;Z_33;X_34;Y_34;Z_34;X_35;Y_35;Z_35;X_36;Y_36;Z_36;X_37;Y_37;Z_37;X_38;Y_38;Z_38;X_39;Y_39;Z_39;X_40;Y_40;Z_40;X_41;Y_41;Z_41;X_42;Y_42;Z_42;X_43;Y_43;Z_43;X_44;Y_44;Z_44;X_45;Y_45;Z_45;X_46;Y_46;Z_46;X_47;Y_47;Z_47;X_48;Y_48;Z_48;X_49;Y_49;Z_49;X_50;Y_50;Z_50;X_51;Y_51;Z_51;X_52;Y_52;Z_52;X_53;Y_53;Z_53;X_54;Y_54;Z_54;X_55;Y_55;Z_55;X_56;Y_56;Z_56;X_57;Y_57;Z_57;X_58;Y_58;Z_58;X_59;Y_59;Z_59;X_60;Y_60;Z_60;X_61;Y_61;Z_61;X_62;Y_62;Z_62;X_63;Y_63;Z_63;X_64;Y_64;Z_64;X_65;Y_65;Z_65;X_66;Y_66;Z_66;X_67;Y_67;Z_67;X_68;Y_68;Z_68;X_69;Y_69;Z_69;X_70;Y_70;Z_70;X_71;Y_71;Z_71;X_72;Y_72;Z_72;X_73;Y_73;Z_73;X_74;Y_74;Z_74;X_75;Y_75;Z_75;X_76;Y_76;Z_76;X_77;Y_77;Z_77;X_78;Y_78;Z_78;X_79;Y_79;Z_79;X_80;Y_80;Z_80;X_81;Y_81;Z_81;X_82;Y_82;Z_82;X_83;Y_83;Z_83;X_84;Y_84;Z_84;X_85;Y_85;Z_85;X_86;Y_86;Z_86;X_87;Y_87;Z_87;X_88;Y_88;Z_88;X_89;Y_89;Z_89;X_90;Y_90;Z_90;X_91;Y_91;Z_91;X_92;Y_92;Z_92;X_93;Y_93;Z_93;X_94;Y_94;Z_94;X_95;Y_95;Z_95;X_96;Y_96;Z_96;X_97;Y_97;Z_97;X_98;Y_98;Z_98;X_99;Y_99;Z_99;X_100;Y_100;Z_100;X_101;Y_101;Z_101;X_102;Y_102;Z_102;X_103;Y_103;Z_103;X_104;Y_104;Z_104;X_105;Y_105;Z_105;X_106;Y_106;Z_106;X_107;Y_107;Z_107;X_108;Y_108;Z_108;X_109;Y_109;Z_109;X_110;Y_110;Z_110;X_111;Y_111;Z_111;X_112;Y_112;Z_112;X_113;Y_113;Z_113;X_114;Y_114;Z_114;X_115;Y_115;Z_115;X_116;Y_116;Z_116;X_117;Y_117;Z_117;X_118;Y_118;Z_118;X_119;Y_119;Z_119;X_120;Y_120;Z_120;X_121;Y_121;Z_121;X_122;Y_122;Z_122;X_123;Y_123;Z_123;X_124;Y_124;Z_124;X_125;Y_125;Z_125;X_126;Y_126;Z_126;X_127;Y_127;Z_127;X_128;Y_128;Z_128;X_129;Y_129;Z_129;X_130;Y_130;Z_130;X_131;Y_131;Z_131;X_132;Y_132;Z_132;X_133;Y_133;Z_133;X_134;Y_134;Z_134;X_135;Y_135;Z_135;X_136;Y_136;Z_136;X_137;Y_137;Z_137;X_138;Y_138;Z_138;X_139;Y_139;Z_139;X_140;Y_140;Z_140;X_141;Y_141;Z_141;X_142;Y_142;Z_142;X_143;Y_143;Z_143;X_144;Y_144;Z_144;X_145;Y_145;Z_145;X_146;Y_146;Z_146;X_147;Y_147;Z_147;X_148;Y_148;Z_148;X_149;Y_149;Z_149;X_150;Y_150;Z_150\n')
	data.close()
	count = 1
	for TheFile in tqdm.tqdm(pdbfilelist):
		data = open(TheFile , 'r')
		seen = set()
		coordinates = list()
		for line in data:
			if line.split()[0] == 'ATOM' and line.split()[2] == 'CA':
				line_lower = line.split()[5].lower()
				if line_lower not in seen:
					seen.add(line_lower)
					line = [char for char in line]
					x = ''.join(line[30:38]).strip()
					y = ''.join(line[38:46]).strip()
					z = ''.join(line[46:54]).strip()
					coordinates.append(x)
					coordinates.append(y)
					coordinates.append(z)
		if len(coordinates) > 450:
			continue
		addition = 450 - len(coordinates)
		for zeros in range(addition):
			coordinates.append('0')
		COline = ';'.join(coordinates)
		data = open('dataCA.csv' , 'a')
		data.write(str(count) + ';' + TheFile.split('.')[0] + ';' + COline + '\n')
		data.close()
		count += 1
	os.chdir(current)
	os.rename(directory + '/dataCA.csv' , 'dataCA.csv')

def DatasetPSO(directory):
	''' Get each residue's phi, psi, and omega angles (uses the PyRosetta library) '''
	''' Generates a the dataPSO.csv with the phi, psi, and omega angles for each amino acid '''
	current = os.getcwd()
	pdbfilelist = os.listdir(directory)
	os.chdir(directory)
	print('\x1b[32m' + "Getting the psi, psi, and omega angles" + '\x1b[0m')
	data = open('dataPSO.csv' , 'a')
	data.write(';PDB_ID;phi_1;psi_1;omg_1;phi_2;psi_2;omg_2;phi_3;psi_3;omg_3;phi_4;psi_4;omg_4;phi_5;psi_5;omg_5;phi_6;psi_6;omg_6;phi_7;psi_7;omg_7;phi_8;psi_8;omg_8;phi_9;psi_9;omg_9;phi_10;psi_10;omg_10;phi_11;psi_11;omg_11;phi_12;psi_12;omg_12;phi_13;psi_13;omg_13;phi_14;psi_14;omg_14;phi_15;psi_15;omg_15;phi_16;psi_16;omg_16;phi_17;psi_17;omg_17;phi_18;psi_18;omg_18;phi_19;psi_19;omg_19;phi_20;psi_20;omg_20;phi_21;psi_21;omg_21;phi_22;psi_22;omg_22;phi_23;psi_23;omg_23;phi_24;psi_24;omg_24;phi_25;psi_25;omg_25;phi_26;psi_26;omg_26;phi_27;psi_27;omg_27;phi_28;psi_28;omg_28;phi_29;psi_29;omg_29;phi_30;psi_30;omg_30;phi_31;psi_31;omg_31;phi_32;psi_32;omg_32;phi_33;psi_33;omg_33;phi_34;psi_34;omg_34;phi_35;psi_35;omg_35;phi_36;psi_36;omg_36;phi_37;psi_37;omg_37;phi_38;psi_38;omg_38;phi_39;psi_39;omg_39;phi_40;psi_40;omg_40;phi_41;psi_41;omg_41;phi_42;psi_42;omg_42;phi_43;psi_43;omg_43;phi_44;psi_44;omg_44;phi_45;psi_45;omg_45;phi_46;psi_46;omg_46;phi_47;psi_47;omg_47;phi_48;psi_48;omg_48;phi_49;psi_49;omg_49;phi_50;psi_50;omg_50;phi_51;psi_51;omg_51;phi_52;psi_52;omg_52;phi_53;psi_53;omg_53;phi_54;psi_54;omg_54;phi_55;psi_55;omg_55;phi_56;psi_56;omg_56;phi_57;psi_57;omg_57;phi_58;psi_58;omg_58;phi_59;psi_59;omg_59;phi_60;psi_60;omg_60;phi_61;psi_61;omg_61;phi_62;psi_62;omg_62;phi_63;psi_63;omg_63;phi_64;psi_64;omg_64;phi_65;psi_65;omg_65;phi_66;psi_66;omg_66;phi_67;psi_67;omg_67;phi_68;psi_68;omg_68;phi_69;psi_69;omg_69;phi_70;psi_70;omg_70;phi_71;psi_71;omg_71;phi_72;psi_72;omg_72;phi_73;psi_73;omg_73;phi_74;psi_74;omg_74;phi_75;psi_75;omg_75;phi_76;psi_76;omg_76;phi_77;psi_77;omg_77;phi_78;psi_78;omg_78;phi_79;psi_79;omg_79;phi_80;psi_80;omg_80;phi_81;psi_81;omg_81;phi_82;psi_82;omg_82;phi_83;psi_83;omg_83;phi_84;psi_84;omg_84;phi_85;psi_85;omg_85;phi_86;psi_86;omg_86;phi_87;psi_87;omg_87;phi_88;psi_88;omg_88;phi_89;psi_89;omg_89;phi_90;psi_90;omg_90;phi_91;psi_91;omg_91;phi_92;psi_92;omg_92;phi_93;psi_93;omg_93;phi_94;psi_94;omg_94;phi_95;psi_95;omg_95;phi_96;psi_96;omg_96;phi_97;psi_97;omg_97;phi_98;psi_98;omg_98;phi_99;psi_99;omg_99;phi_100;psi_100;omg_100;phi_101;psi_101;omg_101;phi_102;psi_102;omg_102;phi_103;psi_103;omg_103;phi_104;psi_104;omg_104;phi_105;psi_105;omg_105;phi_106;psi_106;omg_106;phi_107;psi_107;omg_107;phi_108;psi_108;omg_108;phi_109;psi_109;omg_109;phi_110;psi_110;omg_110;phi_111;psi_111;omg_111;phi_112;psi_112;omg_112;phi_113;psi_113;omg_113;phi_114;psi_114;omg_114;phi_115;psi_115;omg_115;phi_116;psi_116;omg_116;phi_117;psi_117;omg_117;phi_118;psi_118;omg_118;phi_119;psi_119;omg_119;phi_120;psi_120;omg_120;phi_121;psi_121;omg_121;phi_122;psi_122;omg_122;phi_123;psi_123;omg_123;phi_124;psi_124;omg_124;phi_125;psi_125;omg_125;phi_126;psi_126;omg_126;phi_127;psi_127;omg_127;phi_128;psi_128;omg_128;phi_129;psi_129;omg_129;phi_130;psi_130;omg_130;phi_131;psi_131;omg_131;phi_132;psi_132;omg_132;phi_133;psi_133;omg_133;phi_134;psi_134;omg_134;phi_135;psi_135;omg_135;phi_136;psi_136;omg_136;phi_137;psi_137;omg_137;phi_138;psi_138;omg_138;phi_139;psi_139;omg_139;phi_140;psi_140;omg_140;phi_141;psi_141;omg_141;phi_142;psi_142;omg_142;phi_143;psi_143;omg_143;phi_144;psi_144;omg_144;phi_145;psi_145;omg_145;phi_146;psi_146;omg_146;phi_147;psi_147;omg_147;phi_148;psi_148;omg_148;phi_149;psi_149;omg_149;phi_150;psi_150;omg_150\n')
	data.close()
	count = 1
	for TheFile in tqdm.tqdm(pdbfilelist):
		pose = pose_from_pdb(TheFile)
		size = len(pose)
		angles = list()
		for aa in range(size):
			phi = pose.phi(aa + 1)
			psi = pose.psi(aa + 1)
			omg = pose.omega(aa + 1)
			angles.append(str(phi) + ';' + str(psi) + ';' + str(omg))
		Angles = ';'.join(angles)
		if len(angles) >= 150:
			AngLine = Angles
		else:
			addition = 150 - len(angles)
			zeros = list()
			for adds in range(addition):
				zeros.append('0.0;0.0;0.0')
			Zeros = ';'.join(zeros)
			AngLine = Angles + ';' + Zeros
		data = open('dataPSO.csv' , 'a')
		data.write(str(count) + ';' + TheFile + ';' + AngLine + '\n')
		data.close()
		count += 1

def DatasetPS(directory):
	''' Get each residue's phi and psi angles (uses the BioPython library) '''
	''' Generates a the dataPS.csv with the phi and psi angles for each amino acid '''
	current = os.getcwd()
	pdbfilelist = os.listdir(directory)
	os.chdir(directory)
	print('\x1b[32m' + "Getting the psi, psi, and omega angles" + '\x1b[0m')
	data = open('dataPS.csv' , 'a')
	data.write(';PDB_ID;phi_1;psi_1;phi_2;psi_2;phi_3;psi_3;phi_4;psi_4;phi_5;psi_5;phi_6;psi_6;phi_7;psi_7;phi_8;psi_8;phi_9;psi_9;phi_10;psi_10;phi_11;psi_11;phi_12;psi_12;phi_13;psi_13;phi_14;psi_14;phi_15;psi_15;phi_16;psi_16;phi_17;psi_17;phi_18;psi_18;phi_19;psi_19;phi_20;psi_20;phi_21;psi_21;phi_22;psi_22;phi_23;psi_23;phi_24;psi_24;phi_25;psi_25;phi_26;psi_26;phi_27;psi_27;phi_28;psi_28;phi_29;psi_29;phi_30;psi_30;phi_31;psi_31;phi_32;psi_32;phi_33;psi_33;phi_34;psi_34;phi_35;psi_35;phi_36;psi_36;phi_37;psi_37;phi_38;psi_38;phi_39;psi_39;phi_40;psi_40;phi_41;psi_41;phi_42;psi_42;phi_43;psi_43;phi_44;psi_44;phi_45;psi_45;phi_46;psi_46;phi_47;psi_47;phi_48;psi_48;phi_49;psi_49;phi_50;psi_50;phi_51;psi_51;phi_52;psi_52;phi_53;psi_53;phi_54;psi_54;phi_55;psi_55;phi_56;psi_56;phi_57;psi_57;phi_58;psi_58;phi_59;psi_59;phi_60;psi_60;phi_61;psi_61;phi_62;psi_62;phi_63;psi_63;phi_64;psi_64;phi_65;psi_65;phi_66;psi_66;phi_67;psi_67;phi_68;psi_68;phi_69;psi_69;phi_70;psi_70;phi_71;psi_71;phi_72;psi_72;phi_73;psi_73;phi_74;psi_74;phi_75;psi_75;phi_76;psi_76;phi_77;psi_77;phi_78;psi_78;phi_79;psi_79;phi_80;psi_80;phi_81;psi_81;phi_82;psi_82;phi_83;psi_83;phi_84;psi_84;phi_85;psi_85;phi_86;psi_86;phi_87;psi_87;phi_88;psi_88;phi_89;psi_89;phi_90;psi_90;phi_91;psi_91;phi_92;psi_92;phi_93;psi_93;phi_94;psi_94;phi_95;psi_95;phi_96;psi_96;phi_97;psi_97;phi_98;psi_98;phi_99;psi_99;phi_100;psi_100;phi_101;psi_101;phi_102;psi_102;phi_103;psi_103;phi_104;psi_104;phi_105;psi_105;phi_106;psi_106;phi_107;psi_107;phi_108;psi_108;phi_109;psi_109;phi_110;psi_110;phi_111;psi_111;phi_112;psi_112;phi_113;psi_113;phi_114;psi_114;phi_115;psi_115;phi_116;psi_116;phi_117;psi_117;phi_118;psi_118;phi_119;psi_119;phi_120;psi_120;phi_121;psi_121;phi_122;psi_122;phi_123;psi_123;phi_124;psi_124;phi_125;psi_125;phi_126;psi_126;phi_127;psi_127;phi_128;psi_128;phi_129;psi_129;phi_130;psi_130;phi_131;psi_131;phi_132;psi_132;phi_133;psi_133;phi_134;psi_134;phi_135;psi_135;phi_136;psi_136;phi_137;psi_137;phi_138;psi_138;phi_139;psi_139;phi_140;psi_140;phi_141;psi_141;phi_142;psi_142;phi_143;psi_143;phi_144;psi_144;phi_145;psi_145;phi_146;psi_146;phi_147;psi_147;phi_148;psi_148;phi_149;psi_149;phi_150;psi_150\n')
	data.close()
	count = 1
	for TheFile in tqdm.tqdm(pdbfilelist):
		structure = Bio.PDB.PDBParser().get_structure('X' , TheFile)
		dssp = Bio.PDB.DSSP(structure[0] , TheFile , acc_array = 'Wilke')
		angles = list()
		for aa in dssp:
			phi = aa[4]
			psi = aa[5]
			angles.append(str(phi) + ';' + str(psi))
		Angles = ';'.join(angles)
		if len(angles) >= 150:
			AngLine = Angles
		else:
			addition = 150 - len(angles)
			zeros = list()
			for adds in range(addition):
				zeros.append('0.0;0.0')
			Zeros = ';'.join(zeros)
			AngLine = Angles + ';' + Zeros
		data = open('dataPS.csv' , 'a')
		data.write(str(count) + ';' + TheFile + ';' + AngLine + '\n')
		data.close()
		count += 1
		
def DatasetPSOC(directory):
	''' Get each residue's phi, psi, and omega angles as well as CA atom constraints (uses the PyRosetta library) '''
	''' Generates a the dataPSOC.csv with the phi, psi, and omega angles as well as CA atom constraints for each amino acid '''
	current = os.getcwd()
	pdbfilelist = os.listdir(directory)
	os.chdir(directory)
	print('\x1b[32m' + "Getting the psi, psi, and omega angles and CA atom constraints" + '\x1b[0m')
	data = open('dataPSOC.csv' , 'a')
	data.write(';PDB_ID;phi_1;psi_1;omg_1;cst_1;phi_2;psi_2;omg_2;cst_2;phi_3;psi_3;omg_3;cst_3;phi_4;psi_4;omg_4;cst_4;phi_5;psi_5;omg_5;cst_5;phi_6;psi_6;omg_6;cst_6;phi_7;psi_7;omg_7;cst_7;phi_8;psi_8;omg_8;cst_8;phi_9;psi_9;omg_9;cst_9;phi_10;psi_10;omg_10;cst_10;phi_11;psi_11;omg_11;cst_11;phi_12;psi_12;omg_12;cst_12;phi_13;psi_13;omg_13;cst_13;phi_14;psi_14;omg_14;cst_14;phi_15;psi_15;omg_15;cst_15;phi_16;psi_16;omg_16;cst_16;phi_17;psi_17;omg_17;cst_17;phi_18;psi_18;omg_18;cst_18;phi_19;psi_19;omg_19;cst_19;phi_20;psi_20;omg_20;cst_20;phi_21;psi_21;omg_21;cst_21;phi_22;psi_22;omg_22;cst_22;phi_23;psi_23;omg_23;cst_23;phi_24;psi_24;omg_24;cst_24;phi_25;psi_25;omg_25;cst_25;phi_26;psi_26;omg_26;cst_26;phi_27;psi_27;omg_27;cst_27;phi_28;psi_28;omg_28;cst_28;phi_29;psi_29;omg_29;cst_29;phi_30;psi_30;omg_30;cst_30;phi_31;psi_31;omg_31;cst_31;phi_32;psi_32;omg_32;cst_32;phi_33;psi_33;omg_33;cst_33;phi_34;psi_34;omg_34;cst_34;phi_35;psi_35;omg_35;cst_35;phi_36;psi_36;omg_36;cst_36;phi_37;psi_37;omg_37;cst_37;phi_38;psi_38;omg_38;cst_38;phi_39;psi_39;omg_39;cst_39;phi_40;psi_40;omg_40;cst_40;phi_41;psi_41;omg_41;cst_41;phi_42;psi_42;omg_42;cst_42;phi_43;psi_43;omg_43;cst_43;phi_44;psi_44;omg_44;cst_44;phi_45;psi_45;omg_45;cst_45;phi_46;psi_46;omg_46;cst_46;phi_47;psi_47;omg_47;cst_47;phi_48;psi_48;omg_48;cst_48;phi_49;psi_49;omg_49;cst_49;phi_50;psi_50;omg_50;cst_50;phi_51;psi_51;omg_51;cst_51;phi_52;psi_52;omg_52;cst_52;phi_53;psi_53;omg_53;cst_53;phi_54;psi_54;omg_54;cst_54;phi_55;psi_55;omg_55;cst_55;phi_56;psi_56;omg_56;cst_56;phi_57;psi_57;omg_57;cst_57;phi_58;psi_58;omg_58;cst_58;phi_59;psi_59;omg_59;cst_59;phi_60;psi_60;omg_60;cst_60;phi_61;psi_61;omg_61;cst_61;phi_62;psi_62;omg_62;cst_62;phi_63;psi_63;omg_63;cst_63;phi_64;psi_64;omg_64;cst_64;phi_65;psi_65;omg_65;cst_65;phi_66;psi_66;omg_66;cst_66;phi_67;psi_67;omg_67;cst_67;phi_68;psi_68;omg_68;cst_68;phi_69;psi_69;omg_69;cst_69;phi_70;psi_70;omg_70;cst_70;phi_71;psi_71;omg_71;cst_71;phi_72;psi_72;omg_72;cst_72;phi_73;psi_73;omg_73;cst_73;phi_74;psi_74;omg_74;cst_74;phi_75;psi_75;omg_75;cst_75;phi_76;psi_76;omg_76;cst_76;phi_77;psi_77;omg_77;cst_77;phi_78;psi_78;omg_78;cst_78;phi_79;psi_79;omg_79;cst_79;phi_80;psi_80;omg_80;cst_80;phi_81;psi_81;omg_81;cst_81;phi_82;psi_82;omg_82;cst_82;phi_83;psi_83;omg_83;cst_83;phi_84;psi_84;omg_84;cst_84;phi_85;psi_85;omg_85;cst_85;phi_86;psi_86;omg_86;cst_86;phi_87;psi_87;omg_87;cst_87;phi_88;psi_88;omg_88;cst_88;phi_89;psi_89;omg_89;cst_89;phi_90;psi_90;omg_90;cst_90;phi_91;psi_91;omg_91;cst_91;phi_92;psi_92;omg_92;cst_92;phi_93;psi_93;omg_93;cst_93;phi_94;psi_94;omg_94;cst_94;phi_95;psi_95;omg_95;cst_95;phi_96;psi_96;omg_96;cst_96;phi_97;psi_97;omg_97;cst_97;phi_98;psi_98;omg_98;cst_98;phi_99;psi_99;omg_99;cst_99;phi_100;psi_100;omg_100;cst_100;phi_101;psi_101;omg_101;cst_101;phi_102;psi_102;omg_102;cst_102;phi_103;psi_103;omg_103;cst_103;phi_104;psi_104;omg_104;cst_104;phi_105;psi_105;omg_105;cst_105;phi_106;psi_106;omg_106;cst_106;phi_107;psi_107;omg_107;cst_107;phi_108;psi_108;omg_108;cst_108;phi_109;psi_109;omg_109;cst_109;phi_110;psi_110;omg_110;cst_110;phi_111;psi_111;omg_111;cst_111;phi_112;psi_112;omg_112;cst_112;phi_113;psi_113;omg_113;cst_113;phi_114;psi_114;omg_114;cst_114;phi_115;psi_115;omg_115;cst_115;phi_116;psi_116;omg_116;cst_116;phi_117;psi_117;omg_117;cst_117;phi_118;psi_118;omg_118;cst_118;phi_119;psi_119;omg_119;cst_119;phi_120;psi_120;omg_120;cst_120;phi_121;psi_121;omg_121;cst_121;phi_122;psi_122;omg_122;cst_122;phi_123;psi_123;omg_123;cst_123;phi_124;psi_124;omg_124;cst_124;phi_125;psi_125;omg_125;cst_125;phi_126;psi_126;omg_126;cst_126;phi_127;psi_127;omg_127;cst_127;phi_128;psi_128;omg_128;cst_128;phi_129;psi_129;omg_129;cst_129;phi_130;psi_130;omg_130;cst_130;phi_131;psi_131;omg_131;cst_131;phi_132;psi_132;omg_132;cst_132;phi_133;psi_133;omg_133;cst_133;phi_134;psi_134;omg_134;cst_134;phi_135;psi_135;omg_135;cst_135;phi_136;psi_136;omg_136;cst_136;phi_137;psi_137;omg_137;cst_137;phi_138;psi_138;omg_138;cst_138;phi_139;psi_139;omg_139;cst_139;phi_140;psi_140;omg_140;cst_140;phi_141;psi_141;omg_141;cst_141;phi_142;psi_142;omg_142;cst_142;phi_143;psi_143;omg_143;cst_143;phi_144;psi_144;omg_144;cst_144;phi_145;psi_145;omg_145;cst_145;phi_146;psi_146;omg_146;cst_146;phi_147;psi_147;omg_147;cst_147;phi_148;psi_148;omg_148;cst_148;phi_149;psi_149;omg_149;cst_149;phi_150;psi_150;omg_150;cst_150\n')
	data.close()
	count = 1
	for TheFile in tqdm.tqdm(pdbfilelist):
		pyrosetta.toolbox.cleaning.cleanATOM(TheFile)
		TheFile2 = TheFile.split('.')
		pose = pose_from_pdb('{}.clean.pdb'.format(TheFile2[0]))
		size = len(pose)
		phi = list()
		psi = list()
		omg = list()
		cst = list()
		for aa in range(size):
			p = pose.phi(aa + 1)
			#Convert all phi angle values to 0 to 360 (rather than +180 to -180)
			if p < 0:
				p = p + 360
			phi.append(p)
			s = pose.psi(aa + 1)
			#Convert all psi angle values to 0 to 360 (rather than +180 to -180)
			if s < 0:
				s = s + 360
			psi.append(s)
			o = pose.omega(aa + 1)
			#Convert all omega angle values to 0 to 360 (rather than +180 to -180)
			if o < 0:
				o = o + 360
			omg.append(o)
		structure = Bio.PDB.PDBParser().get_structure('X' , TheFile)
		dssp = Bio.PDB.DSSP(structure[0] , TheFile , acc_array = 'Wilke')
		for aa in dssp:
			length = aa[0]
		structure = Bio.PDB.PDBParser(QUIET = True).get_structure('X' , TheFile)
		ppb = Bio.PDB.Polypeptide.PPBuilder()
		Type = ppb.build_peptides(structure , aa_only = False)
		model = Type
		chain = model[0]
		cst.append(0.0)
		for aa in range(1 , length + 1):
			try:
				residue1 = chain[0]
				residue2 = chain[aa]
				atom1 = residue1['CA']
				atom2 = residue2['CA']
				cst.append(atom1 - atom2)
			except:
				pass
		angles = list()
		for P , S , O , C in zip(phi , psi , omg , cst):
			angles.append(str(round(P , 3)) + ';' + str(round(S , 3)) + ';' + str(round(O , 3)) + ';' + str(round(C , 3)))
		Angles = ';'.join(angles)
		if len(angles) >= 150:
			AngLine = Angles
		else:
			addition = 150 - len(angles)
			zeros = list()
			for adds in range(addition):
				zeros.append('0.0;0.0;0.0;0.0')
			Zeros = ';'.join(zeros)
			AngLine = Angles + ';' + Zeros
		TheLine = str(count) + ';' + TheFile + ';' + AngLine + '\n'
		data = open('dataPSOC.csv' , 'a')
		data.write(TheLine)
		data.close()
		count += 1
	os.system('mv dataPSOC.csv {}'.format(current))

def DatasetPSC(directory):
	''' Get each residue's phi and psi angles as well as CA atom constraints (uses the PyRosetta library) '''
	''' Generates a the dataPSC.csv with the phi and psi angles as well as CA atom constraints for each amino acid '''
	current = os.getcwd()
	pdbfilelist = os.listdir(directory)
	os.chdir(directory)
	print('\x1b[32m' + "Getting the psi, psi, and omega angles and CA atom constraints" + '\x1b[0m')
	data = open('dataPSC.csv' , 'a')
	data.write(';PDB_ID;phi_1;psi_1;cst_1;phi_2;psi_2;cst_2;phi_3;psi_3;cst_3;phi_4;psi_4;cst_4;phi_5;psi_5;cst_5;phi_6;psi_6;cst_6;phi_7;psi_7;cst_7;phi_8;psi_8;cst_8;phi_9;psi_9;cst_9;phi_10;psi_10;cst_10;phi_11;psi_11;cst_11;phi_12;psi_12;cst_12;phi_13;psi_13;cst_13;phi_14;psi_14;cst_14;phi_15;psi_15;cst_15;phi_16;psi_16;cst_16;phi_17;psi_17;cst_17;phi_18;psi_18;cst_18;phi_19;psi_19;cst_19;phi_20;psi_20;cst_20;phi_21;psi_21;cst_21;phi_22;psi_22;cst_22;phi_23;psi_23;cst_23;phi_24;psi_24;cst_24;phi_25;psi_25;cst_25;phi_26;psi_26;cst_26;phi_27;psi_27;cst_27;phi_28;psi_28;cst_28;phi_29;psi_29;cst_29;phi_30;psi_30;cst_30;phi_31;psi_31;cst_31;phi_32;psi_32;cst_32;phi_33;psi_33;cst_33;phi_34;psi_34;cst_34;phi_35;psi_35;cst_35;phi_36;psi_36;cst_36;phi_37;psi_37;cst_37;phi_38;psi_38;cst_38;phi_39;psi_39;cst_39;phi_40;psi_40;cst_40;phi_41;psi_41;cst_41;phi_42;psi_42;cst_42;phi_43;psi_43;cst_43;phi_44;psi_44;cst_44;phi_45;psi_45;cst_45;phi_46;psi_46;cst_46;phi_47;psi_47;cst_47;phi_48;psi_48;cst_48;phi_49;psi_49;cst_49;phi_50;psi_50;cst_50;phi_51;psi_51;cst_51;phi_52;psi_52;cst_52;phi_53;psi_53;cst_53;phi_54;psi_54;cst_54;phi_55;psi_55;cst_55;phi_56;psi_56;cst_56;phi_57;psi_57;cst_57;phi_58;psi_58;cst_58;phi_59;psi_59;cst_59;phi_60;psi_60;cst_60;phi_61;psi_61;cst_61;phi_62;psi_62;cst_62;phi_63;psi_63;cst_63;phi_64;psi_64;cst_64;phi_65;psi_65;cst_65;phi_66;psi_66;cst_66;phi_67;psi_67;cst_67;phi_68;psi_68;cst_68;phi_69;psi_69;cst_69;phi_70;psi_70;cst_70;phi_71;psi_71;cst_71;phi_72;psi_72;cst_72;phi_73;psi_73;cst_73;phi_74;psi_74;cst_74;phi_75;psi_75;cst_75;phi_76;psi_76;cst_76;phi_77;psi_77;cst_77;phi_78;psi_78;cst_78;phi_79;psi_79;cst_79;phi_80;psi_80;cst_80;phi_81;psi_81;cst_81;phi_82;psi_82;cst_82;phi_83;psi_83;cst_83;phi_84;psi_84;cst_84;phi_85;psi_85;cst_85;phi_86;psi_86;cst_86;phi_87;psi_87;cst_87;phi_88;psi_88;cst_88;phi_89;psi_89;cst_89;phi_90;psi_90;cst_90;phi_91;psi_91;cst_91;phi_92;psi_92;cst_92;phi_93;psi_93;cst_93;phi_94;psi_94;cst_94;phi_95;psi_95;cst_95;phi_96;psi_96;cst_96;phi_97;psi_97;cst_97;phi_98;psi_98;cst_98;phi_99;psi_99;cst_99;phi_100;psi_100;cst_100;phi_101;psi_101;cst_101;phi_102;psi_102;cst_102;phi_103;psi_103;cst_103;phi_104;psi_104;cst_104;phi_105;psi_105;cst_105;phi_106;psi_106;cst_106;phi_107;psi_107;cst_107;phi_108;psi_108;cst_108;phi_109;psi_109;cst_109;phi_110;psi_110;cst_110;phi_111;psi_111;cst_111;phi_112;psi_112;cst_112;phi_113;psi_113;cst_113;phi_114;psi_114;cst_114;phi_115;psi_115;cst_115;phi_116;psi_116;cst_116;phi_117;psi_117;cst_117;phi_118;psi_118;cst_118;phi_119;psi_119;cst_119;phi_120;psi_120;cst_120;phi_121;psi_121;cst_121;phi_122;psi_122;cst_122;phi_123;psi_123;cst_123;phi_124;psi_124;cst_124;phi_125;psi_125;cst_125;phi_126;psi_126;cst_126;phi_127;psi_127;cst_127;phi_128;psi_128;cst_128;phi_129;psi_129;cst_129;phi_130;psi_130;cst_130;phi_131;psi_131;cst_131;phi_132;psi_132;cst_132;phi_133;psi_133;cst_133;phi_134;psi_134;cst_134;phi_135;psi_135;cst_135;phi_136;psi_136;cst_136;phi_137;psi_137;cst_137;phi_138;psi_138;cst_138;phi_139;psi_139;cst_139;phi_140;psi_140;cst_140;phi_141;psi_141;cst_141;phi_142;psi_142;cst_142;phi_143;psi_143;cst_143;phi_144;psi_144;cst_144;phi_145;psi_145;cst_145;phi_146;psi_146;cst_146;phi_147;psi_147;cst_147;phi_148;psi_148;cst_148;phi_149;psi_149;cst_149;phi_150;psi_150;cst_150\n')
	data.close()
	count = 1
	for TheFile in tqdm.tqdm(pdbfilelist):
		structure = Bio.PDB.PDBParser(QUIET = True).get_structure('X' , TheFile)
		dssp = Bio.PDB.DSSP(structure[0] , TheFile , acc_array = 'Wilke')
		for aa in dssp:
			length = aa[0]
		phi = list()
		psi = list()
		cst = list()
		for aa in dssp:
			#Convert all phi angle values to 0 to 360 (rather than +180 to -180)
			p = aa[4]
			if p < 0:
				p = p + 360
			phi.append(p)
			#Convert all psi angle values to 0 to 360 (rather than +180 to -180)
			s = aa[5]
			if s < 0:
				s = s + 360
			psi.append(s)
		ppb = Bio.PDB.Polypeptide.PPBuilder()
		Type = ppb.build_peptides(structure , aa_only = False)
		model = Type
		chain = model[0]
		cst.append(0.0)
		for aa in range(1 , length + 1):
			try:
				residue1 = chain[0]
				residue2 = chain[aa]
				atom1 = residue1['CA']
				atom2 = residue2['CA']
				cst.append(atom1 - atom2)
			except:
				pass
		angles = list()
		for P , S , C in zip(phi , psi , cst):
			angles.append(str(round(P , 3)) + ';' + str(round(S , 3)) + ';' + str(round(C , 3)))
		Angles = ';'.join(angles)
		if len(angles) >= 150:
			AngLine = Angles
		else:
			addition = 150 - len(angles)
			zeros = list()
			for adds in range(addition):
				zeros.append('0.0;0.0;0.0')
			Zeros = ';'.join(zeros)
			AngLine = Angles + ';' + Zeros
		TheLine = str(count) + ';' + TheFile + ';' + AngLine + '\n'
		data = open('dataPSC.csv' , 'a')
		data.write(TheLine)
		data.close()
		count += 1
	os.system('mv dataPSC.csv {}'.format(current))

#---------------------------------------------------------------------------------------------------------------------------------------
#Protocol to isolate specific types of structures
Database('DATABASE' , 'PDBDatabase')	# 1. Download the PDB database
Extract('PDBDatabase')			# 2. Extract files
NonProtein('PDBDatabase')		# 3. Remove non-protein structures
Size('PDBDatabase' , 80 , 150)		# 4. Remove structures less than or larger than a specified amino acid leangth
Break('PDBDatabase')			# 5. Remove structure with broken chains
Loops('PDBDatabase' , 10)		# 6. Remove structures that have loops that are larger than a spesific length
Renumber('PDBDatabase')			# 7. Renumber structures starting at amino acid 1
#RMSD('PDBDatabase' , 5)		# 8. Measure RMSD of each structure to each structure, remove if RMSD < specified value (CODE IS NOT VERY RELIABLE)
Rg('PDBDatabase' , 15)			# 9. Remove structures that are below a specified Raduis of Gyration value
Sequence('PDBDatabase' , 75)		# 10. Align the sequences of each structure to each structure, remove structures with similar sequences that fall above a user defined percentage

#Protocol to extract specific information from isolated structures
#DatasetR('PDBDatabase')		# 11. Get the secondary structures and distances
#DatasetCA('PDBDatabase')		# 12. Get each residue's CA atom's XYZ coordinates
#DatasetPSO('PDBDatabase')		# 13. Get each residue's phi, psi, and omega angles
#DatasetPS('PDBDatabase')		# 14. Get each residue's phi and psi angles
#DatasetPSOC('PDBDatabase')		# 15. Get each residue's phi, psi, and omega angles as well as CA atom constraints
DatasetPSC('PDBDatabase')		# 16. Get each residue's phi and psi angles as well as CA atom constraints