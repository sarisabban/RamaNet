#!/usr/bin/python3

import os , math , gzip , Bio.PDB

def Database(To , From):
	''' A small script that cleans the PDB database, then isolates the secondary structure and the Phi/Psi torsion angles from each .pdb file '''
	''' Will generate the PDBDatabase directory with all the cleaned .pdb structures inside it, and the Data directory that contains the .csv files for all .pdb files '''
	From = int(smaller)
	To = int(bigger)
	#Collect structures
	os.system('rsync -rlpt -v -z --delete --port=33444 rsync.wwpdb.org::ftp/data/structures/divided/pdb/ ./DATABASE')
	thedatafile = open('data.csv' , 'a')
	thedatafile.write(';1;2;3;4;5;6;7;8;9;10;11;12;13;14;15;16;17;18;19;20;21;22;23;24;25;26;27;28;29;30;31;32;33;34;35;36;37;38;39;40;41;42;43;44;45;46;47;48;49;50;51;52;53;54;55;56;57;58;59;60;61;62;63;64;65;66;67;68;69;70;71;72;73;74;75;76;77;78;79;80;81;82;83;84;85;86;87;88;89;90;91;92;93;94;95;96;97;98;99;100;101;102;103;104;105;106;107;108;109;110;111;112;113;114;115;116;117;118;119;120;121;122;123;124;125;126;127;128;129;130;131;132;133;134;135;136;137;138;139;140;141;142;143;144;145;146;147;148;149;150;Distance_1;Distance_2;Distance_3;Distance_4;Distance_5;Distance_6;Distance_7;Distance_8;Distance_9;Distance_10;\n')
	current = os.getcwd()
	os.mkdir('PDBDatabase')
	filelist = os.listdir('DATABASE')
	for directories in filelist:
		files = os.listdir(current + '/DATABASE/' + directories)
		for afile in files:
			location = (current + '/DATABASE/' + directories + '/' + afile)
			print(location)
			os.rename(location , current + '/PDBDatabase/' + afile)
	os.system('rm -r ./DATABASE')
	#Clean Database
	pdbfilelist = os.listdir('PDBDatabase')
	io = Bio.PDB.PDBIO()
	os.chdir('PDBDatabase')
	for thefile in pdbfilelist:
		try:
			#Open file
			TheFile = current + '/PDBDatabase/' + thefile
			TheName = thefile.split('.')[0].split('pdb')[1].upper()
			#Extract file
			InFile = gzip.open(TheFile, 'rt')
			#Separate chains and save to different files
			structure = Bio.PDB.PDBParser(QUIET=True).get_structure(TheName , InFile)
			count = 0
			for chain in structure.get_chains():
				io.set_structure(chain)
				io.save(structure.get_id() + '_' + chain.get_id() + '.pdb')
			print('\x1b[32m' + '[+] Extracted' + '\t' + thefile.upper() + '\x1b[0m')
			os.remove(TheFile)

		except:
			print('\x1b[31m' + '[-] Failed to Extracted' + '\t' + thefile.upper() + '\x1b[0m')
			os.remove(TheFile)
	os.chdir(current)
	#Remove unwanted structures
	current = os.getcwd()
	pdbfilelist = os.listdir('PDBDatabase')
	ProteinCount = 1
	current = os.getcwd()
	for thefile in pdbfilelist:
		try:
			TheFile = current + '/PDBDatabase/' + thefile
			structure = Bio.PDB.PDBParser(QUIET=True).get_structure('X' , TheFile)
			ppb = Bio.PDB.Polypeptide.PPBuilder()
			Type = ppb.build_peptides(structure , aa_only=True)
			#Delete non-protein files
			if Type == []:
				print('\x1b[31m' + '[-] NOT PROTEIN\t' , thefile + '\x1b[0m')
				os.remove(TheFile)
			else:
				#Delete structures larger than 150 or smaller than 100 amino acids
				parser = Bio.PDB.PDBParser()
				structure = parser.get_structure('X' , TheFile)
				model = structure[0]
				dssp = Bio.PDB.DSSP(model , TheFile , acc_array='Wilke')
				for aa in dssp:
					length = aa[0]
				if length >= To or length <= From:
					print('\x1b[31m' + '[-] WRONG SIZE\t' , thefile + '\x1b[0m')
					os.remove(TheFile)
				else:
					#Delete structures with none-continuous chains by tracing the chain and measuring all the peptide bonds (aprox = 1.3 angstroms), if the distance between the C and N atoms is larger than 1.3 then there is a chain break
					structure = Bio.PDB.PDBParser(QUIET=True).get_structure('X' , TheFile)
					ppb = Bio.PDB.Polypeptide.PPBuilder()
					Type = ppb.build_peptides(structure , aa_only=True)
					model = Type
					chain = model[0]
					count = 0
					ChainBreak = None
					for bond in range(length):
						try:
							residue1 = chain[count]
							residue2 = chain[count + 1]
							atom1 = residue1['C']
							atom2 = residue2['N']
							distance = atom1-atom2
							count += 1
							if distance > 1.4:
								ChainBreak = 'Break'
							else:
								pass
						except:
							pass
					if ChainBreak == 'Break':
						print('\x1b[31m' + '[-] CHAIN BREAK\t' , thefile + '\x1b[0m')
						os.remove(TheFile)
					else:
						#Get secondary structures
						parser = Bio.PDB.PDBParser()
						structure = parser.get_structure('X' , TheFile)
						model = structure[0]
						dssp = Bio.PDB.DSSP(model , TheFile , acc_array='Wilke')
						SS = list()
						for res in dssp:
							ss = res[2]
							if ss == '-' or ss == 'T' or ss == 'S':		#Loop (DSSP code is - or T or S)
								SS.append('L')
							elif ss == 'G' or ss == 'H' or ss == 'I':	#Helix (DSSP code is G or H or I)
								SS.append('H')
							elif ss == 'B' or ss == 'E':			#Sheet (DSSP code is B or E)
								SS.append('S')
						#Delete floppy structures, with loops as their dominant secondary structure
						loop = SS.count('L')
						helix = SS.count('H')
						strand = SS.count('S')
						if loop >= helix + strand:
							print('\x1b[31m' + '[-] FLOPPY\t' , thefile + '\x1b[0m')
							os.remove(TheFile)
						else:
							#Calculate Rg
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
										mass.append(32.065)
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
							if rg <= 15:
								print('\x1b[31m' + '[-] HIGH Rg\t' , thefile + '\x1b[0m')
								os.remove(TheFile)
							else:
								#Renumber residues
								pdb = open(TheFile , 'r')
								PDB = open(TheFile + 'X' , 'w')
								count = 0
								num = 0
								AA2 = None
								for line in pdb:
									count += 1														#Sequencially number atoms
									AA1 = line[23:27]													#Sequencially number residues
									if not AA1 == AA2:
										num += 1			
									final_line = line[:7] + '{:4d}'.format(count) + line[11:17] + line[17:21] + 'A' + '{:4d}'.format(num) + line[26:]	#Update each line to have its atoms and residues sequencially labeled, as well as being in chain A
									AA2 = AA1
									PDB.write(final_line)													#Write to new file called motif.pdb
								PDB.close()
								os.remove(TheFile)
								os.rename(TheFile + 'X' , TheFile)
								#Get torsion angles
								count += 1
								Tor = list()
								for model in Bio.PDB.PDBParser().get_structure('X' , TheFile):
									for chain in model:
										polypeptides = Bio.PDB.PPBuilder().build_peptides(chain)
										for poly_index , poly in enumerate(polypeptides):
											phi_psi = poly.get_phi_psi_list()
											for res_index , residue in enumerate(poly):
												#Phi angles
												if phi_psi[res_index][0] is None:
													phi = 0
												else:
													angle = phi_psi[res_index][0] * 180 / math.pi
													while angle > 180:
														angle = angle - 360
													while angle < -180:
														angle = angle + 360
													phi = angle
												#Psi angles
												if phi_psi[res_index][1] is None:
													psi = 0
												else:
													angle = phi_psi[res_index][1] * 180 / math.pi
													while angle > 180:
														angle = angle - 360
													while angle < -180:
														angle = angle + 360
													psi = angle
												Tor.append((phi , psi))
								#Distances
								structure = Bio.PDB.PDBParser(QUIET=True).get_structure('X' , TheFile)
								ppb = Bio.PDB.Polypeptide.PPBuilder()
								Type = ppb.build_peptides(structure , aa_only=False)
								model = Type
								chain = model[0]
								length = int(str(Type[0]).split()[2].split('=')[1].split('>')[0])
								distances = list()
								if length >= To or length <= From:
									print('\x1b[31m' + '[-] WRONG SIZE\t' , thefile + '\x1b[0m')
									os.remove(TheFile)
								else:
									for key , value in {0:1 , 9:11 , 19:21 , 29:31 , 39:41 , 49:51 , 59:61 , 69:71 , 79:81 , 89:91}.items():
										residue1 = chain[key]
										residue2 = chain[length - value]
										atom1 = residue1['CA']
										atom2 = residue2['CA']
										distance = atom1-atom2
										distances.append(distance)
									#Put info together
									ss = list()
									for val in SS:
										if val == 'L':
											ss.append('1')
										elif val == 'H':
											ss.append('2')
										elif val == 'S':
											ss.append('3')
									SecondaryStructures = ';' + ';'.join(ss)					#Secondary Structures L = 1, H = 2, S = 3 printed horisantally
									phiang = list()
									psiang = list()
									for val in Tor:
										phiang.append(val[0])
										psiang.append(val[1])
									PHIAngles = ';' + ';'.join(map(str, phiang))					#PHI angles printed horisantally
									PSIAngles = ';' + ';'.join(map(str, psiang))					#PSI angles printed horisantally
									Distances = ';' + ';'.join(map(str , distances))
									#Fill in remaining positions with 0 until position 150
									add = 150 - len(SS)
									fill = list()
									for zeros in range(add):
										fill.append('0')
									filling = ';' + ';'.join(fill)
									#Write to file
									line = str(ProteinCount) + SecondaryStructures + filling + Distances + '\n'	#The PHI and PSI angels are not being used because we cannot insert the angels as a feature during Machine Learning prediction, to use add this to the line variable: PHIAngles + filling + PSIAngles + filling
									thedatafile.write(line)
									ProteinCount += 1
									print('\x1b[32m' + '[+] GOOD\t' , thefile + '\x1b[0m')
		except:
			print('\x1b[31m' + '[-] Script Crashed' + '\t' + thefile.upper() + '\x1b[0m')
	thedatafile.close()
	os.system('rm -r PDBDatabase')

Database(100 , 150)
