#!/usr/bin/python3

import os , re , time , datetime , random , requests , urllib.request , bs4 , Bio.PDB
from pyrosetta import *
from pyrosetta.toolbox import *
init()
#--------------------------------------------------------------------------------------------------------------------------------------
#Functions

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
	score1_original_before_relax = scorefxn(pose)										#Measure score before relaxing
	scorefxn = get_fa_scorefxn()
	relax = pyrosetta.rosetta.protocols.relax.FastRelax()
	relax.set_scorefxn(scorefxn)
	relax.apply(pose)
	score2_original_after_relax = scorefxn(pose)										#Measure score after relaxing
	#B - FastDesign protocol												#Uses Generic Monte Carlo with PackStat as a filter to direct FastDesign towards an optimally packed structure core
	chain = pose.pdb_info().chain(1)											#Identify chain
	layers = [2 , 1 , 0]													#Layer Identity from SASA Surface = [0] , Boundary = [1] , Core = [2]
	for identity in layers:													#Loop through each layer
		#1 - Setup the PackStat filter
		filters = rosetta.protocols.simple_filters.PackStatFilter()
		#2 - Identify The Layers
		sasa = SASA(pose)												#Re-calculate SASA every time because amino acid position can change from one layer to another during the design phase, therefore make sure to design the layer not the amino acid
		layer = sasa[identity]												#Changes every iteration to start with Core (sasa[2]) then Boundary (sasa[1]) then Surface (sasa[0])
		#3 - Generate the resfile											#Will generate a new Resfile for each layer (which is why it is deleted at the end of the loop)
		Resfile = open('Resfile.resfile' , 'w')
		Resfile.write('NATAA\n')
		Resfile.write('start\n')
		for line in layer:
			Resfile.write(str(line) + ' ' + chain + ' ALLAA\n')
		Resfile.close()
		#4 - Setup the FastDesign mover
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
		#5 - Setup and apply the generic Monte Carlo mover
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
	#C - Relax pose
	relax.apply(pose)
	#D - Output result
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

def Draw(SecondaryStructureString , DistancesList):
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
	#Generate constraints file
	ConstFile = open('constraints.cst' , 'w')
	firstAA = 0
	secndAA = length
	for distance in DistancesList:
		if firstAA == 0:
			line = 'AtomPair CA ' + '1' + ' CA ' + str(secndAA) + ' GAUSSIANFUNC ' + str(distance) + ' 2.0\n'
			ConstFile.write(line)
			firstAA += 10
			secndAA -= 10
		else:
			line = 'AtomPair CA ' + str(firstAA) + ' CA ' + str(secndAA) + ' GAUSSIANFUNC ' + str(distance) + ' 2.0\n'
			ConstFile.write(line)
			firstAA += 10
			secndAA -= 10
	ConstFile.close()
	#Fold topology
	constraints = pyrosetta.rosetta.protocols.simple_moves.ConstraintSetMover()
	constraints.constraint_file('constraints.cst')
	constraints.add_constraints(True)
	constraints.apply(pose)



#	X = pyrosetta.rosetta.protocols.relax.AtomCoordinateCstMover()
#	X.set_type('AtomPair CA 1 CA 124 GAUSSIANFUNC 1.0 1.0')
#	X.apply(pose)


	scorefxn = get_fa_scorefxn()
	relax = pyrosetta.rosetta.protocols.relax.FastRelax()
	relax.set_scorefxn(scorefxn)
	relax.constrain_relax_to_start_coords(True)
	relax.constrain_coords(True)
	relax.apply(pose)



	pose.dump_pdb('DeNovo.pdb')
	os.remove('constraints.cst')




















def GenSecStruct():
	''' Generates a random structure's secondary structures '''
	''' Returns a string that has each amino acid's secondary structure '''
	#The secondary structure ratios
	Ratio_H_to_S_Number = [2 , 5]###########################################MUST STUDY THE RATIO OF HELIX TO STRAND TO LOOP ---> WHAT IS THE PATTERN???
	Ratio_H_to_S_Size = [10 , 10]###########################################MUST STUDY THE RATIO OF HELIX TO STRAND TO LOOP ---> WHAT IS THE PATTERN???
	#Generate Random protein size
	size = random.randint(100 , 150)
	helix_number = Ratio_H_to_S_Number[0]
	helix_size = Ratio_H_to_S_Size[0]
	strand_number = Ratio_H_to_S_Number[1]
	strand_size = Ratio_H_to_S_Size[1]
	#Generate helices
	Helix = list()
	for numb in range(helix_number):
		Hsize = list()
		for H in range(helix_size):
			Hsize.append('H')
		Ahelix = ''.join(Hsize)
		Helix.append(Ahelix)
	#Generate strands
	Strand = list()
	for numb in range(strand_number):
		Ssize = list()
		for S in range(strand_size):
			Ssize.append('S')
		Astrand = ''.join(Ssize)
		Strand.append(Astrand)
	#Generate loops
	while True:
		order = Helix + Strand
		random.shuffle(order)
		HandS_size = (helix_size * helix_number) + (strand_size * strand_number)
		loop_number = len(order) - 1
		loop_size = size - HandS_size - 2
		Ls = list()
		for numb in range(loop_size):
			Ls.append('L')
		for chunk in range(loop_number):
			Ls.append('.')
		random.shuffle(Ls)
		Ls = ''.join(Ls)
		Loop = Ls.split('.')
		decision = None
		for check in Loop:
			if len(check) < 3:
				decision = 'Bad'
				break
			else:
				decision = 'Good'
		if decision == 'Good':
			break
		else:
			continue
	#Put together
	protein = ['L']
	count = 0
	while True:
		try:
			protein.append(order[count])
			protein.append(Loop[count])
			count += 1
		except:
			break
	protein.append('L')
	protein = ''.join(protein)

	return(protein)










def ML(CSV_FILE):
	''' Takes the data.csv and learns the distances pattern given the secondary structure of each amino acid, this is to allow it to generate distances in an effort to fold a topology into a logical protein structure '''
	''' Returns a list of the distances between spesific parts of the protein that can be used as constrains when drawing and folding the DeNovo protein's topology '''
	pass
#--------------------------------------------------------------------------------------------------------------------------------------
'''
SASA(pose)
Design(pose)
Fragments(pose)
SS = GenSecStruct()
dist = ML('Data.csv')
Draw(SS , dist)
'''
#--------------------------------------------------------------------------------------------------------------------------------------
SS = 'LLLHHHHHHHHHLLLLLLLLLHHHHHHHHHHHLLLLLLLLLLSSSSSLLLHHHHHHHHHLSSSLLLLLLLLSSSLLLLSSSSSSSSLLLLSLLSLLSSSSSSSSLSSSSSSLLLLSSSSSSSSL'
dist = [46.2503 , 40.3013 , 25.9238 , 11.769 , 11.1069 , 16.1582 , 12.93 , 18.1924 , 14.2343 , 16.1879]
Draw(SS , dist)
