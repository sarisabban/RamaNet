#!/usr/bin/python3

import sys , os , re , time , datetime , subprocess , random , requests , itertools , urllib.request , bs4 , math , gzip , Bio.PDB.Polypeptide , Bio.PDB.PDBParser , Bio.PDB, Bio.PDB.Vector
from pyrosetta import *
from pyrosetta.toolbox import *
init()

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
		ss = open('structure.blueprint' , 'w')
		for residues in range(size):
			for x in info:
				if residues == x[0]:
					for y in range(x[1]):
						ss.write('L' + '\n')				#Loop insert
			ss.write('H' + '\n')							#Helix insert
		ss.close()

def Database(filename):
	''' Isolates the secondary structure and the Phi/Psi torsion angles from a .pdb file '''
	''' Generates a .dat file for that .pdb file '''
	#Get secondary structures
	parser = Bio.PDB.PDBParser()
	structure = parser.get_structure('structure' , filename)
	model = structure[0]
	dssp = Bio.PDB.DSSP(model , filename , acc_array='Wilke')
	SS = list()
	for res in dssp:
		ss = res[2]
		if ss == '-' or ss == 'T' or ss == 'S':		#Loop (DSSP code is - or T or S)
			SS.append('L')
		elif ss == 'G' or ss == 'H' or ss == 'I':	#Helix (DSSP code is G or H or I)
			SS.append('H')
		elif ss == 'B' or ss == 'E':			#Sheet (DSSP code is B or E)
			SS.append('S')
	#Get torsion angles
	Tor = list()
	for model in Bio.PDB.PDBParser().get_structure('structure' , filename):
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
	#Put together
	name = filename.split('.')
	thefile = open(name[0] + '.data' , 'w')
	for item in zip(SS , Tor):
		line = (str(item[0]) + '\t' + str(item[1][0]) + '\t' + str(item[1][1]) + '\n')
		thefile.write(line)
	thefile.close()







def PDBclean(smaller , bigger):
	''' A small script that cleans the PDB database '''
	''' Will generate the directory PDBDatabase with all the cleaned .pdb structures inside it '''
	From = int(smaller)
	To = int(bigger)

	#Collect Structures
	os.system('rsync -rlpt -v -z --delete --port=33444 rsync.wwpdb.org::ftp/data/structures/divided/pdb/ ./DATABASE')
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
	for thefile in pdbfilelist:
		TheFile = current + '/PDBDatabase/' + thefile
		structure = Bio.PDB.PDBParser(QUIET=True).get_structure('X' , TheFile)
		ppb = Bio.PDB.Polypeptide.PPBuilder()
		Type = ppb.build_peptides(structure)
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
			else:
	print('[+] GOOD\t' , thefile)
















































def DeNovo(filename):
	''' Draws the torsion angles to generate a PDB file '''
	''' Generates the Draw.pdb file '''
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




def Learn():
	pass
#--------------------------------------------------------------------------------------------------------------------------------------
Database('start.pdb')
Draw('start.data')
pose = pose_from_pdb('DeNovo.pdb')
Design(pose)
Fragments(pose)
