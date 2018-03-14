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

def DatasetPS(directory):
	''' Get each residue's phi and psi angles '''
	''' Generates a the dataPS.csv with the XYZ coordinates of the CA atom for each amino acid '''
	current = os.getcwd()
	pdbfilelist = os.listdir(directory)
	os.chdir(directory)
	print('\x1b[32m' + "Getting the psi and psi angles" + '\x1b[0m')
	data = open('dataCA.csv' , 'a')
	data.write(';PDB_ID;phi_1;psi_1;phi_2;psi_2;phi_3;psi_3;phi_4;psi_4;phi_5;psi_5;phi_6;psi_6;phi_7;psi_7;phi_8;psi_8;phi_9;psi_9;phi_10;psi_10;phi_11;psi_11;phi_12;psi_12;phi_13;psi_13;phi_14;psi_14;phi_15;psi_15;phi_16;psi_16;phi_17;psi_17;phi_18;psi_18;phi_19;psi_19;phi_20;psi_20;phi_21;psi_21;phi_22;psi_22;phi_23;psi_23;phi_24;psi_24;phi_25;psi_25;phi_26;psi_26;phi_27;psi_27;phi_28;psi_28;phi_29;psi_29;phi_30;psi_30;phi_31;psi_31;phi_32;psi_32;phi_33;psi_33;phi_34;psi_34;phi_35;psi_35;phi_36;psi_36;phi_37;psi_37;phi_38;psi_38;phi_39;psi_39;phi_40;psi_40;phi_41;psi_41;phi_42;psi_42;phi_43;psi_43;phi_44;psi_44;phi_45;psi_45;phi_46;psi_46;phi_47;psi_47;phi_48;psi_48;phi_49;psi_49;phi_50;psi_50;phi_51;psi_51;phi_52;psi_52;phi_53;psi_53;phi_54;psi_54;phi_55;psi_55;phi_56;psi_56;phi_57;psi_57;phi_58;psi_58;phi_59;psi_59;phi_60;psi_60;phi_61;psi_61;phi_62;psi_62;phi_63;psi_63;phi_64;psi_64;phi_65;psi_65;phi_66;psi_66;phi_67;psi_67;phi_68;psi_68;phi_69;psi_69;phi_70;psi_70;phi_71;psi_71;phi_72;psi_72;phi_73;psi_73;phi_74;psi_74;phi_75;psi_75;phi_76;psi_76;phi_77;psi_77;phi_78;psi_78;phi_79;psi_79;phi_80;psi_80;phi_81;psi_81;phi_82;psi_82;phi_83;psi_83;phi_84;psi_84;phi_85;psi_85;phi_86;psi_86;phi_87;psi_87;phi_88;psi_88;phi_89;psi_89;phi_90;psi_90;phi_91;psi_91;phi_92;psi_92;phi_93;psi_93;phi_94;psi_94;phi_95;psi_95;phi_96;psi_96;phi_97;psi_97;phi_98;psi_98;phi_99;psi_99;phi_100;psi_100;phi_101;psi_101;phi_102;psi_102;phi_103;psi_103;phi_104;psi_104;phi_105;psi_105;phi_106;psi_106;phi_107;psi_107;phi_108;psi_108;phi_109;psi_109;phi_110;psi_110;phi_111;psi_111;phi_112;psi_112;phi_113;psi_113;phi_114;psi_114;phi_115;psi_115;phi_116;psi_116;phi_117;psi_117;phi_118;psi_118;phi_119;psi_119;phi_120;psi_120;phi_121;psi_121;phi_122;psi_122;phi_123;psi_123;phi_124;psi_124;phi_125;psi_125;phi_126;psi_126;phi_127;psi_127;phi_128;psi_128;phi_129;psi_129;phi_130;psi_130;phi_131;psi_131;phi_132;psi_132;phi_133;psi_133;phi_134;psi_134;phi_135;psi_135;phi_136;psi_136;phi_137;psi_137;phi_138;psi_138;phi_139;psi_139;phi_140;psi_140;phi_141;psi_141;phi_142;psi_142;phi_143;psi_143;phi_144;psi_144;phi_145;psi_145;phi_146;psi_146;phi_147;psi_147;phi_148;psi_148;phi_149;psi_149;phi_150;psi_150\n')
	data.close()
	count = 1
	for TheFile in tqdm.tqdm(pdbfilelist):





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
DatasetR('PDBDatabase')			# 11. Get the secondary structures and distances
DatasetCA('PDBDatabase')		# 12. Get each residue's CA atom's XYZ coordinates
DatasetPS('PDBDatabase')		# 13. Get each residue's phi and psi angles
