#!/usr/bin/python3

import os , re , time , datetime , random , requests , urllib.request , bs4 , math , gzip , Bio.PDB
from pyrosetta import *
from pyrosetta.toolbox import *
init()
#--------------------------------------------------------------------------------------------------------------------------------------
#Functions

def Relax(pose):
	''' Relaxes a structure '''
	''' Updates the original pose with the relaxed pose '''
	scorefxn = get_fa_scorefxn()
	relax = pyrosetta.rosetta.protocols.relax.FastRelax()
	relax.set_scorefxn(scorefxn)
	relax.apply(pose)

def SASA(pose):
	''' Calculates the different layers (Surface, Boundery, Core) of a structure according its SASA (solvent-accessible surface area) '''
	''' Returns three lists Surface amino acids = [0] , Boundery amino acids = [1] , Core amino acids = [2] '''
	#Temporary generate a .pdb file of the pose to isolate the layers since it is not yet possible to do that using a Rosetta pose, this temporary .pdb file will be deleted after the layers are found
	pose.dump_pdb('ToDesign.pdb')
	#Standard script to setup biopython's DSSP to calculate SASA using Wilke constants
	parser = Bio.PDB.PDBParser()
	structure = parser.get_structure('X' , 'ToDesign.pdb')
	model = structure[0]
	dssp = Bio.PDB.DSSP(model , 'ToDesign.pdb' , acc_array='Wilke')
	#Loop to get SASA for each amino acid
	lis = list()
	count = 0
	for x in dssp:
		if x[1]=='A' : sasa=129*(x[3])
		elif x[1]=='V' : sasa=174*(x[3])
		elif x[1]=='I' : sasa=197*(x[3])
		elif x[1]=='L' : sasa=201*(x[3])
		elif x[1]=='M' : sasa=224*(x[3])
		elif x[1]=='P' : sasa=159*(x[3])
		elif x[1]=='Y' : sasa=263*(x[3])
		elif x[1]=='F' : sasa=240*(x[3])
		elif x[1]=='W' : sasa=285*(x[3])
		elif x[1]=='R' : sasa=274*(x[3])
		elif x[1]=='C' : sasa=167*(x[3])
		elif x[1]=='N' : sasa=195*(x[3])
		elif x[1]=='Q' : sasa=225*(x[3])
		elif x[1]=='E' : sasa=223*(x[3])
		elif x[1]=='G' : sasa=104*(x[3])
		elif x[1]=='H' : sasa=224*(x[3])
		elif x[1]=='K' : sasa=236*(x[3])
		elif x[1]=='S' : sasa=155*(x[3])
		elif x[1]=='T' : sasa=172*(x[3])
		elif x[1]=='D' : sasa=193*(x[3])
		lis.append((x[2] , sasa))
	#Label each amino acid depending on its SASA position according to the parameters highlighted in the paper by (Koga et.al., 2012 - PMID: 23135467). The parameters are as follows:
	#Surface:	Helix or Sheet: SASA => 60		Loop: SASA => 40
	#Boundry:	Helix or Sheet: 15 < SASA < 60		Loop: 25 < SASA < 40
	#Core:		Helix or Sheet: SASA =< 15		Loop: SASA =< 25	
	surface = list()
	boundery = list()
	core = list()
	count = 0
	for x , y in lis:
		count = count + 1
		if y <= 25 and (x == '-' or x == 'T' or x == 'S'):		#Loop (DSSP code is - or T or S)
			core.append(count)
		elif 25 < y < 40 and (x == '-' or x == 'T' or x == 'S'):	#Loop (DSSP code is - or T or S)
			boundery.append(count)
		elif y >= 40 and (x == '-' or x == 'T' or x == 'S'):		#Loop (DSSP code is - or T or S)
			surface.append(count)
		elif y <= 15 and (x == 'G' or x == 'H' or x == 'I'):		#Helix (DSSP code is G or H or I)
			core.append(count)
		elif 15 < y < 60 and (x == 'G' or x == 'H' or x == 'I'):	#Helix (DSSP code is G or H or I)
			boundery.append(count)
		elif y >= 60 and (x == 'G' or x == 'H' or x == 'I'):		#Helix (DSSP code is G or H or I)
			surface.append(count)
		elif y <= 15 and (x == 'B' or x == 'E'):			#Sheet (DSSP code is B or E)
			core.append(count)
		elif 15 < y < 60 and (x == 'B' or x == 'E'):			#Sheet (DSSP code is B or E)
			boundery.append(count)
		elif y >= 60 and (x == 'B' or x == 'E'):			#Sheet (DSSP code is B or E)
			surface.append(count)	
	os.remove('ToDesign.pdb')						#Keep working directory clean
	return(surface , boundery , core)													#Return values [0] = Motif_From [1] = Motif_To

def Design(pose):
	''' Applies FastDesign to change the whole structure's amino acids (one layer at a time as well as designing towards an optimally packed core) while maintaining the same backbone. Should be faster than the Whole method and results in a better final structure than the Layer method '''
	''' Generates the Designed.pdb file '''
	#A - Relax original structure
	scorefxn = get_fa_scorefxn()												#Call the score function
	score1_original_before_relax = scorefxn(pose)										#Measure score before relaxing
	Relax(pose)														#Relax structure
	score2_original_after_relax = scorefxn(pose)										#Measure score after relaxing
	#B - FastDesign Protocol												#Uses Generic Monte Carlo with PackStat as a filter to direct FastDesign towards an optimally packed structure core
	chain = pose.pdb_info().chain(1)											#Identify chain
	layers = [2 , 1 , 0]													#Layer Identity from SASA Surface = [0] , Boundary = [1] , Core = [2]
	for identity in layers:													#Loop through each layer
		#1 - Setup The PackStat Filter
		filters = rosetta.protocols.simple_filters.PackStatFilter()
		#2 - Identify The Layers
		sasa = SASA(pose)												#Re-calculate SASA every time because amino acid position can change from one layer to another during the design phase, therefore make sure to design the layer not the amino acid
		layer = sasa[identity]												#Changes every iteration to start with Core (sasa[2]) then Boundary (sasa[1]) then Surface (sasa[0])
		#3 - Generate The Resfile											#Will generate a new Resfile for each layer (which is why it is deleted at the end of the loop)
		Resfile = open('Resfile.resfile' , 'w')
		Resfile.write('NATAA\n')
		Resfile.write('start\n')
		for line in layer:
			Resfile.write(str(line) + ' ' + chain + ' ALLAA\n')
		Resfile.close()
		#4 - Setup The FastDesign Mover
		task = pyrosetta.rosetta.core.pack.task.TaskFactory()								#Setup the TaskFactory
		read = pyrosetta.rosetta.core.pack.task.operation.ReadResfile('Resfile.resfile')				#Call the generated Resfile
		task.push_back(read)												#Add the Resfile to the TaskFactory
		movemap = MoveMap()												#Setup the MoveMap
		movemap.set_bb(False)												#Do not change the phi and psi BackBone angles
		movemap.set_chi(True)												#Change the chi Side Chain angle
		mover = pyrosetta.rosetta.protocols.denovo_design.movers.FastDesign()						#Call the FastDesign Mover
		mover.set_task_factory(task)											#Add the TaskFactory to it
		mover.set_movemap(movemap)											#Add the MoveMap to it
		mover.set_scorefxn(scorefxn)											#Add the Score Function to it
		#5 - Setup and Apply The Generic Monte Carlo Mover
		MC = pyrosetta.rosetta.protocols.simple_moves.GenericMonteCarloMover()						#Call Monter Carlo Class
		MC.set_mover(mover)												#Load The Mover
		MC.set_scorefxn(scorefxn)											#Set score function
		MC.set_maxtrials(10)												#Set number of monte carlo loops
		MC.set_temperature(1)												#Set temperature
		MC.set_preapply(True)												#To apply Boltzmann accept/reject to all applications of the mover (always use False)
		MC.set_drift(True)												#Make current pose = next iteration pose
		MC.set_sampletype('high')											#Move monte carlo to higher filter score
		MC.add_filter(filters , False , 1.0 , 'high' , True)								#Add a filter (Filter Type , Adaptive , Temperature , Sample Type , Rank By)
		MC.apply(pose)													#Apply Move
		os.remove('Resfile.resfile')											#To keep working directory clean, and to make sure each Resfile has the info for each layer only and they do not get mixed and appended together in one Resfile
	#C - Relax Pose
	Relax(pose)														#Relax structure
	#D - Output Result
	score3_of_design_after_relax = scorefxn(pose)										#Measure score of designed pose
	pose.dump_pdb('structure.pdb')												#Export final pose into a .pdb structure file
	print('---------------------------------------------------------')
	print('Original Structure Score:' , '\t' , score1_original_before_relax)
	print('Relaxed Original Score:' , '\t' , score2_original_after_relax)
	print('Relaxed Design Score:' , '\t\t' , score3_of_design_after_relax)

def Fragments(pose):
	''' Submits the pose to the Robetta server (http://www.robetta.org) for fragment generation that are used for the Abinitio folding simulation. Then measures the RMSD for each fragment at each position and chooses the lowest RMSD. Then averages out the lowest RMSDs. Then plots the lowest RMSD fragment for each positon '''
	''' Generates the 3-mer file, the 9-mer file, the PsiPred file, the RMSD vs Position PDF plot with the averaged fragment RMSD printed in the plot '''
	#Make the 3-mer and 9-mer fragment files and the PSIPRED file using the Robetta server
	sequence = pose.sequence()
	#Post
	web = requests.get('http://www.robetta.org/fragmentsubmit.jsp')
	payload = {
		'UserName':'ac.research',
		'Email':'',
		'Notes':'structure',
		'Sequence':sequence,
		'Fasta':'',
		'Code':'',
		'ChemicalShifts':'',
		'NoeConstraints':'',
		'DipolarConstraints':'',
		'type':'submit'
	}
	session = requests.session()
	response = session.post('http://www.robetta.org/fragmentsubmit.jsp', data=payload , files=dict(foo='bar'))		
	for line in response:
		line = line.decode()
		if re.search('<a href="(fragmentqueue.jsp\?id=[0-9].*)">' , line):
			JobID = re.findall('<a href="(fragmentqueue.jsp\?id=[0-9].*)">' , line)
	JobURL = 'http://www.robetta.org/' + JobID[0]
	#Check
	ID = JobID[0].split('=')
	print('Job ID: ' + str(ID[1]))
	while True:
		Job = urllib.request.urlopen(JobURL)
		jobdata = bs4.BeautifulSoup(Job , 'lxml')
		status = jobdata.find('td', string='Status: ').find_next().text
		if status == 'Complete':
			print(datetime.datetime.now().strftime('%d %B %Y @ %H:%M') , 'Status:' , status)
			break
		else:
			print(datetime.datetime.now().strftime('%d %B %Y @ %H:%M') , 'Status:' , status)
			time.sleep(1800)
			continue
	#Download
	sequence = pose.sequence()
	fasta = open('structure.fasta' , 'w')
	fasta.write(sequence)
	fasta.close()
	os.system('wget http://www.robetta.org/downloads/fragments/' + str(ID[1])  + '/aat000_03_05.200_v1_3')
	os.system('wget http://www.robetta.org/downloads/fragments/' + str(ID[1])  + '/aat000_09_05.200_v1_3')
	os.system('wget http://www.robetta.org/downloads/fragments/' + str(ID[1])  + '/t000_.psipred_ss2')
	os.rename('aat000_03_05.200_v1_3' , 'frags.200.3mers')
	os.rename('aat000_09_05.200_v1_3' , 'frags.200.9mers')
	os.rename('t000_.psipred_ss2' , 'pre.psipred.ss2')
	#Calculate the best fragment's RMSD at each position
	frag = open('frags.200.9mers' , 'r')
	rmsd = open('temp.dat' , 'w')
	for line in frag:
		if line.lstrip().startswith('position:'):
			line = line.split()
			size = line[1]
	frag.close()
	count = 0
	for x in range (int(size)):
		count +=1
		#Get the pose and make a copy of it to apply changes to
		pose_copy = Pose()
		pose_copy.assign(pose)
		#Setup frame list
		frames = pyrosetta.rosetta.core.fragment.FrameList()
		#Setup the 9-mer fragment (9-mer is better than 3-mer for this analysis)
		fragset = pyrosetta.rosetta.core.fragment.ConstantLengthFragSet(9)
		fragset.read_fragment_file('frags.200.9mers')
		fragset.frames(count , frames)
		#Setup the MoveMap
		movemap = MoveMap()
		movemap.set_bb(True)
		#Setup and apply the fragment inserting mover
		for frame in frames:
			for frag_num in range( 1 , frame.nr_frags() + 1 ):
				frame.apply(movemap , frag_num , pose_copy)
				#Measure the RMSD difference between the original pose and the new changed pose (the copy)
				RMSD = rosetta.core.scoring.CA_rmsd(pose , pose_copy)
				print(RMSD , '\t' , count)
				rmsd.write(str(RMSD) + '\t' + str(count) + '\n')
				#Reset the copy pose to original pose
				pose_copy.assign(pose)
	rmsd.close()
	#Analyse the RMSD file to get the lowest RMSD for each position
	data = open('RMSDvsPosition.dat' , 'w')
	lowest = {} 									#Mapping group number -> lowest value found
	for line in open('temp.dat'):
		parts = line.split()
		if len(parts) != 2:							#Only lines with two items on it
			continue
		first = float(parts[0])
		second = int(parts[1])
		if first == 0: 								#Skip line with 0.0 RMSD (this is an error from the 9-mer fragment file). I don't know why it happens
			continue
		if second not in lowest:
			lowest[second] = first
		else:
			if first < lowest[second]:
				lowest[second] = first
	for position, rmsd in lowest.items():
		#print(str(rmsd) + '\t' + str(position))
		data.write(str(position) + '\t' + str(rmsd) + '\n')
	data.close()
	#Calculate the average RMSD of the fragments
	data = open('RMSDvsPosition.dat' , 'r')
	value = 0
	for line in data:
		line = line.split()
		RMSD = float(line[1])
		value = value + RMSD
		count = int(line[0])
	Average_RMSD = round(value / count , 2)
	#Plot the results
	gnuplot = open('gnuplot_sets' , 'w')
	gnuplot.write("reset\nset terminal postscript\nset output './plot_frag.pdf'\nset encoding iso_8859_1\nset term post eps enh color\nset xlabel 'Position'\nset ylabel 'RMSD (\\305)'\nset yrange [0:]\nset xrange [0:]\nset xtics auto\nset xtics rotate\nset grid front\nunset grid\nset title 'Fragment Quality'\nset key off\nset boxwidth 0.5\nset style fill solid\nset label 'Average RMSD = " + str(Average_RMSD) + "' at graph 0.01 , graph 0.95 tc lt 7 font 'curior 12'\nplot 'RMSDvsPosition.dat' with boxes\nexit")
	gnuplot.close()
	os.system('gnuplot < gnuplot_sets')
	os.remove('gnuplot_sets')
	os.remove('temp.dat')
	return(Average_RMSD)

def Database(smaller , bigger):
	''' A small script that cleans the PDB database, then isolates the secondary structure and the Phi/Psi torsion angles from each .pdb file '''
	''' Will generate the PDBDatabase directory with all the cleaned .pdb structures inside it, and the Data directory that contains the .csv files for all .pdb files '''
	From = int(smaller)
	To = int(bigger)

	#Collect Structures
	os.system('wget -rA .ent.gz ftp://ftp.rcsb.org/pub/pdb/data/structures/divided/pdb/ -P DATABASE')
	current = os.getcwd()
	os.mkdir('PDBDatabase')
	filelist = os.listdir('DATABASE/ftp.rcsb.org/pub/pdb/data/structures/divided/pdb')
	for directories in filelist:
		files = os.listdir(current + '/DATABASE/ftp.rcsb.org/pub/pdb/data/structures/divided/pdb/' + directories)
		for afile in files:
			location = (current + '/DATABASE/ftp.rcsb.org/pub/pdb/data/structures/divided/pdb/' + directories + '/' + afile)
			print(location)
			os.rename(location , current + '/PDBDatabase/' + afile)
	os.system('rm -r ./DATABASE')
	#Separate Chains
	pdbfilelist = os.listdir('PDBDatabase')
	for thefile in pdbfilelist:
		#Open File
		TheFile = current + '/PDBDatabase/' + thefile
		TheName = TheFile.split('.')
		#Extract Each Chain and Save as Different Files
		InFile = gzip.open(TheFile, 'rb')
		for line in InFile:
			line = line.decode()
			if line.startswith('ATOM') or line.startswith('ANISOU'):
				chain = line[21]
				Name = TheName[0].split('pdb')
				output = open(Name[1] + '_' + chain + '.pdb' , 'a')
				output.write(line)
				output.close()
		os.remove(TheFile)
	#Remove Unwanted Structures
	pdbfilelist = os.listdir('PDBDatabase')
	count = 0
	for thefile in pdbfilelist:
		TheFile = current + '/PDBDatabase/' + thefile
		structure = Bio.PDB.PDBParser(QUIET=True).get_structure('X' , TheFile)
		ppb = Bio.PDB.Polypeptide.PPBuilder()
		Type = ppb.build_peptides(structure , aa_only=False)
		#Delete Non-Protein Files
		if Type == []:
			print('[-] NOT PROTEIN\t' , thefile)
			os.remove(TheFile)
		else:
			#Delete Structures Larger Than 150 or Smaller Than 100 Amino Acids
			length = int(str(Type[0]).split()[2].split('=')[1].split('>')[0])
			if length > To or length < From:
				print('[-] WRONG SIZE\t' , thefile)
				os.remove(TheFile)
			elif:
				#Delete Structures With None-continuous Chains By Tracing The Chain And Measuring All The Peptide Bonds (aprox = 1.3 angstroms), If The Distance Between The C and N Atoms is Larger Than 1.3 Then There Is A Chain Break
				structure = Bio.PDB.PDBParser(QUIET=True).get_structure('X' , TheFile)
				ppb = Bio.PDB.Polypeptide.PPBuilder()
				Type = ppb.build_peptides(structure , aa_only=False)
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
					os.remove(TheFile)
				else:
					pass
			else:
				#Get Secondary Structures
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
				#Delete Floppy Structures, With Loops as Their Dominant Secondary Structure
				loop = SS.count('L')
				helix = SS.count('H')
				strand = SS.count('S')
				if loop >= helix + strand:
					os.remove(TheFile)
				else:
					#Renumber Residues
					pdb = open(TheFile , 'r')
					PDB = open('X' + TheFile , 'w')
					count = 0
					num = 0
					AA2 = None
					for line in pdb:
						count += 1					#Sequencially number atoms
						AA1 = line[23:27]				#Sequencially number residues
						if not AA1 == AA2:
							num += 1			
						final_line = line[:7] + '{:4d}'.format(count) + line[11:17] + line[17:21] + 'A' + '{:4d}'.format(num) + line[26:]	#Update each line of the motif to have its atoms and residues sequencially labeled, as well as being in chain A
						AA2 = AA1
						PDB.write(final_line)				#Write to new file called motif.pdb
					PDB.close()
					os.remove(TheFile)
					os.rename('X' + TheFile , TheFile)
					#Get Torsion Angles
					count += 1
					Tor = list()
					for model in Bio.PDB.PDBParser().get_structure('X' , TheFile):
						for chain in model:
							polypeptides = Bio.PDB.PPBuilder().build_peptides(chain)
							for poly_index , poly in enumerate(polypeptides):
								phi_psi = poly.get_phi_psi_list()
								for res_index , residue in enumerate(poly):
									#Phi Angles
									if phi_psi[res_index][0] is None:
										phi = 0
									else:
										angle = phi_psi[res_index][0] * 180 / math.pi
										while angle > 180:
											angle = angle - 360
										while angle < -180:
											angle = angle + 360
										phi = angle
									#Psi Angles
									if phi_psi[res_index][1] is None:
										psi = 0
									else:
										angle = phi_psi[res_index][1] * 180 / math.pi
										while angle > 180:
											angle = angle - 360
										while angle < -180:
											angle = angle + 360
									Tor.append((phi , psi))
									#Distances
									structure = Bio.PDB.PDBParser(QUIET=True).get_structure('X' , TheFile)
									ppb = Bio.PDB.Polypeptide.PPBuilder()
									Type = ppb.build_peptides(structure , aa_only=False)
									model = Type
									chain = model[0]
									distances = list()
									for key , value in {0:1 , 9:11 , 19:21 , 29:31 , 39:41 , 49:51 , 59:61 , 69:71}.items(): #79:81 , 89:91
										residue1 = chain[key]
										residue2 = chain[length - value]
										atom1 = residue1['CA']
										atom2 = residue2['CA']
										distance = atom1-atom2
										distances.append(distance)
				#Put Together
				name = thefile.split('.')
				thefile = open('data' + '.csv' , 'a')
				thefile.write(';1;2;3;4;5;6;7;8;9;10;11;12;13;14;15;16;17;18;19;20;21;22;23;24;25;26;27;28;29;30;31;32;33;34;35;36;37;38;39;40;41;42;43;44;45;46;47;48;49;50;51;52;53;54;55;56;57;58;59;60;61;62;63;64;65;66;67;68;69;70;71;72;73;74;75;76;77;78;79;80;81;82;83;84;85;86;87;88;89;90;91;92;93;94;95;96;97;98;99;100;101;102;103;104;105;106;107;108;109;110;111;112;113;114;115;116;117;118;119;120;121;122;123;124;125;126;127;128;129;130;131;132;133;134;135;136;137;138;139;140;141;142;143;144;145;146;147;148;149;150;Distance_1;Distance_2;Distance_3;Distance_4;Distance_5;Distance_6;Distance_7;Distance_8;Distance_9;Distance_10;\n')
				#thefile.write(';1;2;3;4;5;6;7;8;9;10;11;12;13;14;15;16;17;18;19;20;21;22;23;24;25;26;27;28;29;30;31;32;33;34;35;36;37;38;39;40;41;42;43;44;45;46;47;48;49;50;51;52;53;54;55;56;57;58;59;60;61;62;63;64;65;66;67;68;69;70;71;72;73;74;75;76;77;78;79;80;81;82;83;84;85;86;87;88;89;90;91;92;93;94;95;96;97;98;99;100;101;102;103;104;105;106;107;108;109;110;111;112;113;114;115;116;117;118;119;120;121;122;123;124;125;126;127;128;129;130;131;132;133;134;135;136;137;138;139;140;141;142;143;144;145;146;147;148;149;150;PHI_1;PHI_2;PHI_3;PHI_4;PHI_5;PHI_6;PHI_7;PHI_8;PHI_9;PHI_10;PHI_11;PHI_12;PHI_13;PHI_14;PHI_15;PHI_16;PHI_17;PHI_18;PHI_19;PHI_20;PHI_21;PHI_22;PHI_23;PHI_24;PHI_25;PHI_26;PHI_27;PHI_28;PHI_29;PHI_30;PHI_31;PHI_32;PHI_33;PHI_34;PHI_35;PHI_36;PHI_37;PHI_38;PHI_39;PHI_40;PHI_41;PHI_42;PHI_43;PHI_44;PHI_45;PHI_46;PHI_47;PHI_48;PHI_49;PHI_50;PHI_51;PHI_52;PHI_53;PHI_54;PHI_55;PHI_56;PHI_57;PHI_58;PHI_59;PHI_60;PHI_61;PHI_62;PHI_63;PHI_64;PHI_65;PHI_66;PHI_67;PHI_68;PHI_69;PHI_70;PHI_71;PHI_72;PHI_73;PHI_74;PHI_75;PHI_76;PHI_77;PHI_78;PHI_79;PHI_80;PHI_81;PHI_82;PHI_83;PHI_84;PHI_85;PHI_86;PHI_87;PHI_88;PHI_89;PHI_90;PHI_91;PHI_92;PHI_93;PHI_94;PHI_95;PHI_96;PHI_97;PHI_98;PHI_99;PHI_100;PHI_101;PHI_102;PHI_103;PHI_104;PHI_105;PHI_106;PHI_107;PHI_108;PHI_109;PHI_110;PHI_111;PHI_112;PHI_113;PHI_114;PHI_115;PHI_116;PHI_117;PHI_118;PHI_119;PHI_120;PHI_121;PHI_122;PHI_123;PHI_124;PHI_125;PHI_126;PHI_127;PHI_128;PHI_129;PHI_130;PHI_131;PHI_132;PHI_133;PHI_134;PHI_135;PHI_136;PHI_137;PHI_138;PHI_139;PHI_140;PHI_141;PHI_142;PHI_143;PHI_144;PHI_145;PHI_146;PHI_147;PHI_148;PHI_149;PHI_150;PHI_1;PSI_2;PSI_3;PSI_4;PSI_5;PSI_6;PSI_7;PSI_8;PSI_9;PSI_10;PSI_11;PSI_12;PSI_13;PSI_14;PSI_15;PSI_16;PSI_17;PSI_18;PSI_19;PSI_20;PSI_21;PSI_22;PSI_23;PSI_24;PSI_25;PSI_26;PSI_27;PSI_28;PSI_29;PSI_30;PSI_31;PSI_32;PSI_33;PSI_34;PSI_35;PSI_36;PSI_37;PSI_38;PSI_39;PSI_40;PSI_41;PSI_42;PSI_43;PSI_44;PSI_45;PSI_46;PSI_47;PSI_48;PSI_49;PSI_50;PSI_51;PSI_52;PSI_53;PSI_54;PSI_55;PSI_56;PSI_57;PSI_58;PSI_59;PSI_60;PSI_61;PSI_62;PSI_63;PSI_64;PSI_65;PSI_66;PSI_67;PSI_68;PSI_69;PSI_70;PSI_71;PSI_72;PSI_73;PSI_74;PSI_75;PSI_76;PSI_77;PSI_78;PSI_79;PSI_80;PSI_81;PSI_82;PSI_83;PSI_84;PSI_85;PSI_86;PSI_87;PSI_88;PSI_89;PSI_90;PSI_91;PSI_92;PSI_93;PSI_94;PSI_95;PSI_96;PSI_97;PSI_98;PSI_99;PSI_100;PSI_101;PSI_102;PSI_103;PSI_104;PSI_105;PSI_106;PSI_107;PSI_108;PSI_109;PSI_110;PSI_111;PSI_112;PSI_113;PSI_114;PSI_115;PSI_116;PSI_117;PSI_118;PSI_119;PSI_120;PSI_121;PSI_122;PSI_123;PSI_124;PSI_125;PSI_126;PSI_127;PSI_128;PSI_129;PSI_130;PSI_131;PSI_132;PSI_133;PSI_134;PSI_135;PSI_136;PSI_137;PSI_138;PSI_139;PSI_140;PSI_141;PSI_142;PSI_143;PSI_144;PSI_145;PSI_146;PSI_147;PSI_148;PSI_149;PSI_150;Distance_1;Distance_2;Distance_3;Distance_4;Distance_5;Distance_6;Distance_7;Distance_8;Distance_9;Distance_10;\n')
				ss = list()
				for val in SS:
					if val == 'L':
						ss.append('1')
					elif val == 'H':
						ss.append('2')
					elif val == 'S':
						ss.append('3')
				SecondaryStructures = ';' + ';'.join(ss)			#Secondary Structures L = 1, H = 2, S = 3 printed horisantally
				phiang = list()
				psiang = list()
				for val in Tor:
					phiang.append(val[0])
					psiang.append(val[1])
				PHIAngles = ';' + ';'.join(map(str, phiang))			#PHI angles printed horisantally
				PSIAngles = ';' + ';'.join(map(str, psiang))			#PSI angles printed horisantally
				Distances = ';' + ';'.join(map(str , distances))
				#Fill in Remaining Positions With 0 Until Position 150
				add = 150 - len(SS)
				fill = list()
				for zeros in range(add):
					fill.append('0')
				filling = ';' + ';'.join(fill)
				line = str(count) + SecondaryStructures + filling + Distances	#The PHI and PSI angels are not being used because we cannot insert the angels as a feature during Machine Learning prediction, to use add this to the line variable: PHIAngles + filling + PSIAngles + filling
				thefile.write(line)
				thefile.close()
				print('[+] GOOD\t' , name[0])
	os.system('mv ' + name[0] + '.csv .')
	os.remove('PDBDatabase')

def Draw(filename):
	''' Draws the torsion angles to generate a .pdb file '''
	''' Generates the DeNovo.pdb file '''
	#Length of structure
	thefile = open(filename , 'r')
	for resi in enumerate(thefile):
		size = resi[0] + 1
	Val = str()
	for itr in range(size):
		itr = 'V'
		Val = Val + itr
	pose = pose_from_sequence(Val)
	#Apply torsion angles
	thefile = open(filename , 'r')
	for line in enumerate(thefile):
		angles = line[1].split()
		phi = angles[1]
		psi = angles[2]
		pose.set_phi(int(line[0] + 1) , float(phi))
		pose.set_psi(int(line[0] + 1) , float(psi))
	Relax(pose)
	pose.dump_pdb('DeNovo.pdb')

def BluePrint():
		''' Generates a random blueprint file '''
		''' Generates the blueprint file '''
		#Generate blueprint file
		size = random.randint(120 , 130)						#Random protein size
		#Construct the loops
		info = list()
		for number in range(random.randint(1 , 4)):					#Randomly choose weather to add 3, 4, or 5 different loops
			Loop = random.randint(0 , 1)						#Randomly choose weather to add a 3 residue loop or a 4 residue loop
			if Loop == 0:
				position = random.randint(1 , size)				#Randomly choose where these loops are added
				info.append((position , 3))
			else:
				position = random.randint(1 , size)
				info.append((position , 4))
		#Generate the blueprint file
		ss = open('blueprint' , 'w')
		for residues in range(size):
			for x in info:
				if residues == x[0]:
					for y in range(x[1]):
						ss.write('L' + '\n')				#Loop insert
			ss.write('H' + '\n')							#Helix insert
		ss.close()








def ML(directory):
	pass
#--------------------------------------------------------------------------------------------------------------------------------------
'''
Relax(pose)
SASA(pose)
Design(pose)
Fragments(pose)
Database(smaller , bigger)
Draw(filename)
BluePrint()
ML(Data)
'''
#--------------------------------------------------------------------------------------------------------------------------------------
Database(100 , 150)
