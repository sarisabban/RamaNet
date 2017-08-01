#!/usr/bin/python3
                               #Modules to download---> #      #         #
import sys , os , re , time , datetime , subprocess , numpy , Bio.PDB , bs4 , random , requests , urllib.request

#Import PyRosetta, its Tools, and its Database
from pyrosetta import *
from pyrosetta.toolbox import *
from pyrosetta.rosetta.protocols.relax import *
from pyrosetta.rosetta.protocols.moves import *
from pyrosetta.rosetta.protocols.simple_moves import *
from pyrosetta.rosetta.protocols.denovo_design.movers import *
from pyrosetta.rosetta.protocols.toolbox.pose_metric_calculators import *
from pyrosetta.rosetta.core.fragment import *
from pyrosetta.rosetta.core.pack.task import *
from pyrosetta.rosetta.core.scoring.packstat import *
from pyrosetta.rosetta.core.pack.task.operation import *
from pyrosetta.rosetta.core.pack.task import TaskFactory
init()

#RosettaRelax
def Relax(pose):
	''' Relaxes a structure '''
	''' Updates the original pose with the relaxed pose '''
	scorefxn = get_fa_scorefxn()
	relax = FastRelax()
	relax.set_scorefxn(scorefxn)
	relax.apply(pose)

#SASA
def SASA(pose):
	''' Calculates the different layers (Surface, Boundery, Core) of a structure according its SASA (solvent-accessible surface area) '''
	''' Returns three lists Surface amino acids = [0] , Boundery amino acids = [1] , Core amino acids = [2] '''
	#Temporary generate a .pdb file of the pose to isolate the layers since it is not yet possible to do that using a Rosetta pose, this temporary .pdb file will be deleted after the layers are found
	pose.dump_pdb('ToDesign.pdb')
	#Standard script to setup biopython's DSSP to calculate SASA using Wilke constants
	p = Bio.PDB.PDBParser()
	structure = p.get_structure('X' , 'ToDesign.pdb')
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
	return(surface , boundery , core)

#RosettaDesign
class Design():
	#Design Whole Structure All At Once
	def Whole(pose):
		''' Applies RosettaDesign to change the whole structure's amino acids (the whole structure all at once) while maintaining the same backbone '''
		''' Just updates the pose with the new structure '''
		#1 - Relax original structure
		scorefxn = get_fa_scorefxn()							#Call the score function
		score1_original_before_relax = scorefxn(pose)					#Measure score before relaxing
		Relax(pose)									#Relax structure
		score2_original_after_relax = scorefxn(pose)					#Measure score after relaxing
		#2 - Preform RosettaDesign for whole structure
		for inter in range(3):
			task_pack = standard_packer_task(pose)
			pack_mover = PackRotamersMover(scorefxn , task_pack)
			pack_mover.apply(pose)
			#3 - Relax pose
			Relax(pose)
		score3_of_design_after_relax = scorefxn(pose)					#Measure score of designed pose
		pose.dump_pdb('structure.pdb')							#Export final pose into a .pdb structure file
		print(score1_original_before_relax)
		print(score2_original_after_relax)
		print(score3_of_design_after_relax)
	#Design The Structure One Layer At A Time Moving Towards A Tightly Packed Core With Every Loop
	def Pack(pose):
		''' Applies FastDesign to change the whole structure's amino acids (one layer at a time as well as designing towards an optimally packed core) while maintaining the same backbone. Should be faster than the Whole method and results in a better final structure than the Layer method '''
		''' Generates the Designed.pdb file '''
		#A - Relax original structure
		scorefxn = get_fa_scorefxn()							#Call the score function
		score1_original_before_relax = scorefxn(pose)					#Measure score before relaxing
		Relax(pose)									#Relax structure
		score2_original_after_relax = scorefxn(pose)					#Measure score after relaxing
		#B - FastDesign Protocol							#Uses Generic Monte Carlo with PackStat as a filter to direct FastDesign towards an optimally packed structure core
		layers = [2 , 1 , 0]								#Layer Identity from SASA Surface = [0] , Boundary = [1] , Core = [2]
		for identity in layers:								#Loop through each layer
			#1 - Setup The PackStat Filter
			filters = rosetta.protocols.simple_filters.PackStatFilter()
			#2 - Identify The Layers
			sasa = SASA(pose)							#Re-calculate SASA every time because amino acid position can change from one layer to another during the design phase, therefore make sure to design the layer not the amino acid
			layer = sasa[identity]							#Changes every iteration to start with Core (sasa[2]) then Boundary (sasa[1]) then Surface (sasa[0])
			#3 - Generate The Resfile						#Will generate a new Resfile for each layer (which is why it is deleted at the end of the loop)
			Resfile = open('Resfile.resfile' , 'w')
			Resfile.write('NATAA\n')
			Resfile.write('start\n')
			for line in layer:
				Resfile.write(str(line) + ' A ALLAA\n')
			Resfile.close()
			#4 - Setup The FastDesign Mover
			task = TaskFactory()							#Setup the TaskFactory
			read = ReadResfile('Resfile.resfile')					#Call the generated Resfile
			task.push_back(read)							#Add the Resfile to the TaskFactory
			movemap = MoveMap()							#Setup the MoveMap
			movemap.set_bb(False)							#Do not change the phi and psi BackBone angles
			movemap.set_chi(True)							#Change the chi Side Chain angle
			mover = FastDesign()							#Call the FastDesign Mover
			mover.set_task_factory(task)						#Add the TaskFactory to it
			mover.set_movemap(movemap)						#Add the MoveMap to it
			mover.set_scorefxn(scorefxn)						#Add the Score Function to it
			#5 - Setup and Apply The Generic Monte Carlo Mover
			MC = GenericMonteCarloMover()						#Call Monter Carlo Class
			MC.set_mover(mover)							#Load The Mover
			MC.set_scorefxn(scorefxn)						#Set score function
			MC.set_maxtrials(1)							#Set number of monte carlo loops
			MC.set_temperature(1)							#Set temperature
			MC.set_preapply(True)							#To apply Boltzmann accept/reject to all applications of the mover (always use False)
			MC.set_drift(True)							#Make current pose = next iteration pose
			MC.set_sampletype('high')						#Move monte carlo to higher filter score
			#MC.recover_low(True)							#True - at the end of application, the pose is set to the lowest (or highest if sample_type="high") scoring pose
			#MC.stopping_condition()						#Stops before trials are done if a filter evaluates to true
			MC.add_filter(filters , False , 1.0 , 'high' , True)			#Add a filter (Filter Type , Adaptive , Temperature , Sample Type , Rank By)
			#MC.task_factory(task) #Causes an infinite loop				#Include a Task Factory
			#MC.boltzmann(pose) #For some reason hates a relaxed pose		#Evaulates a pose based on the scores/filters + temperatures
			MC.apply(pose)								#Apply Move
			os.remove('Resfile.resfile')						#To keep working directory clean, and to make sure each Resfile has the info for each layer only and they do not get mixed and appended together in one Resfile
		#C - Relax Pose
		Relax(pose)									#Relax structure
		#D - Output Result
		score3_of_design_after_relax = scorefxn(pose)					#Measure score of designed pose
		pose.dump_pdb('structure.pdb')							#Export final pose into a .pdb structure file
		print('---------------------------------------------------------')
		print('Original Structure Score:' , '\t' , score1_original_before_relax)
		print('Relaxed Original Score:' , '\t' , score2_original_after_relax)
		print('Relaxed Design Score:' , '\t\t' , score3_of_design_after_relax)

#Fragment Generation and Identification
class Fragment():
	#Make The 3-mer and 9-mer Fragment Files
	def Make(pose):
		''' Submits the pose to the Robetta Server (http://www.robetta.org) for fragment generation that are used for the Abinitio folding simulation '''
		''' Generates the 3-mer file, the 9-mer file, and the PsiPred file '''
		sequence = pose.sequence()
		#1 - Post
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
		#2 - Check
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
		#3 - Download
		sequence = pose.sequence()
		fasta = Resfile = open('structure.fasta' , 'w')
		fasta.write(sequence)
		fasta.close()
		print('wget http://www.robetta.org/downloads/fragments/' + str(ID[1])  + '/aat000_03_05.200_v1_3')
		print('wget http://www.robetta.org/downloads/fragments/' + str(ID[1])  + '/aat000_09_05.200_v1_3')
		print('wget http://www.robetta.org/downloads/fragments/' + str(ID[1])  + '/t000_.psipred_ss2')
		os.system('wget http://www.robetta.org/downloads/fragments/' + str(ID[1])  + '/aat000_03_05.200_v1_3')
		os.system('wget http://www.robetta.org/downloads/fragments/' + str(ID[1])  + '/aat000_09_05.200_v1_3')
		os.system('wget http://www.robetta.org/downloads/fragments/' + str(ID[1])  + '/t000_.psipred_ss2')
	#Calculate The Best Fragment's RMSD At Each Position And Plot The Result
	def RMSD(pose , Fragment_File):
		''' Measures the RMSD for each fragment at each position and plots the lowest RMSD fragment for each positon '''
		''' Generates an RMSD vs Position PDF plot '''
		frag = open(Fragment_File , 'r')
		rmsd = open('temp.dat' , 'w')
		for line in frag:
			if line.startswith(' position:'):
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
			frames = FrameList()
			#Setup the 9-mer fragment (9-mer is better than 3-mer for this analysis)
			fragset = ConstantLengthFragSet(9)
			fragset.read_fragment_file(Fragment_File)
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
		gnuplot = open('gnuplot_sets' , 'w')
		gnuplot.write("""set terminal postscript
set output './plot_frag.pdf'
set encoding iso_8859_1
set term post eps enh color
set xlabel 'Position'
set ylabel 'RMSD (\\305)'
set yrange [0:]
set xrange [0:]
set xtics 1
set xtics rotate 
set title 'Fragment Quality'
set key off
set boxwidth 0.5
set style fill solid
plot 'RMSDvsPosition.dat' with boxes
exit""")
		gnuplot.close()
		os.system('gnuplot < gnuplot_sets')
		os.remove('gnuplot_sets')
		os.remove('temp.dat')

#Denovo Design
def DeNovo(number_of_output):
	''' Preforms De Novo Design on a protein's structure using the BluePrintBDR Mover. Generates only structures with helices (no sheet) '''
	''' Generates user defined number of DeNovo_#.pdb files each with a different structure '''
	#1 - Generate a temporary dummy .pdb so the BluePrintBDR mover can work on
	temp = open('temp.pdb' , 'w')
	temp.write('ATOM      1  N   ASP C   1      33.210  65.401  53.583  1.00 55.66           N  \nATOM      2  CA  ASP C   1      33.590  64.217  54.411  1.00 55.66           C  \nATOM      3  C   ASP C   1      33.574  62.950  53.542  1.00 52.88           C  \nATOM      4  O   ASP C   1      34.516  62.724  52.780  1.00 50.94           O  \nATOM      5  CB  ASP C   1      32.656  64.090  55.624  1.00 58.39           C  ')
	temp.close()
	pose = pose_from_pdb('temp.pdb')
	RgValue = 9999999999
	PoseFinal = Pose()
	for nterm in range(number_of_output):							#Number of different output structures
		#2 - Generate blueprint file
		size = random.randint(120 , 130)						#Random protein size
		#3 - Construct the loops
		info = list()
		for number in range(random.randint(1 , 4)):					#Randomly choose weather to add 3, 4, or 5 different loops
			Loop = random.randint(0 , 1)						#Randomly choose weather to add a 3 residue loop or a 4 residue loop
			if Loop == 0:
				position = random.randint(1 , size)				#Randomly choose where these loops are added
				info.append((position , 3))
			else:
				position = random.randint(1 , size)
				info.append((position , 4))
		#4 - Generate the blueprint file
		blueprint = open('blueprint' , 'w')
		for residues in range(size):
			for x in info:
				if residues == x[0]:
					for y in range(x[1]):
						blueprint.write('0 V ' + 'L' + 'X R\n')		#Loop insert
			blueprint.write('0 V ' + 'H' + 'X R\n')					#Helix insert
		blueprint.close()
		#5 - Run the BluePrintBDR mover
		mover = pyrosetta.rosetta.protocols.fldsgn.BluePrintBDR()
		mover.num_fragpick(200)
		mover.use_fullmer(True)
		mover.use_abego_bias(True)
		mover.use_sequence_bias(False)
		mover.max_linear_chainbreak(0.07)
		mover.ss_from_blueprint(True)
		mover.dump_pdb_when_fail('')
		mover.set_constraints_NtoC(-1.0)
		mover.set_blueprint('blueprint')
		mover.apply(pose)
		os.remove('blueprint')
		#Calculate Radius of Gyration (Rg) and choose lowest Rg score
		Rg = ScoreFunction()
		Rg.set_weight(pyrosetta.rosetta.core.scoring.rg , 1)
		Value = Rg(pose)
		print(Value , '<---------------------------------------')##################
		if Value <= RgValue:
			RgValue = Value
			PoseFinal.assign(pose)
		else:
			continue
	PoseFinal.dump_pdb('DeNovo.pdb')
	os.remove('temp.pdb')
#--------------------------------------------------------------------------------------------------------------------------------------
for structures in range(2):
	directory = os.getcwd()
	folder = 'structure_' + str(structures + 1)
	os.mkdir(folder)
	os.chdir(folder)
	#Start Protocol
	DeNovo(100)
	pose = pose_from_pdb('DeNovo.pdb')
	Design.Pack(pose)
	Fragment.Make(pose)
	Fragment.RMSD(pose , 'aat000_09_05.200_v1_3')
	#End Protocol
	os.chdir(directory)
