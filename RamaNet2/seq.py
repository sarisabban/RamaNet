import os
import tqdm
import math
import Bio
import Bio.PDB
from pyrosetta import *
from pyrosetta.toolbox import *
init('-out:level 0 -ignore_zero_occupancy false')
#assert Bio.__version__ == '1.72', 'Biopython must be version 1.72'

class Dataset():
	def Database(self, TempDIR, FinalDIR):
		'''
		Downloads the entire PDB database from https://www.wwpdb.org/
		moves all files into one directory, then uncompresses all the files
		Generates a directory which contains all .PDB structure files
		'''
		os.system('rsync -rlpt -v -z --delete --port=33444 rsync.wwpdb.org::ftp/data/structures/divided/pdb/ ./{}'.format(TempDIR))
		os.mkdir(FinalDIR)
		filelist = os.listdir(TempDIR)
		print('\x1b[32m[+] Download complete\x1b[0m')
		print('\x1b[32m[+] Moving files\x1b[0m')
		for directories in tqdm.tqdm(filelist):
			files = os.listdir('{}/{}'.format(TempDIR, directories))
			for afile in files:
				location = ('{}/{}/{}'.format(TempDIR, directories, afile))
				os.rename(location, '{}/{}'.format(FinalDIR, afile))
		os.system('rm -r ./{}'.format(TempDIR))
	def Extract(self, directory):
		'''
		Extracts all the .ent.gz files and separate all chains and save them into
		seperate .pdb files. Replaces each .ent.gz file with the .pdb file of each
		chain
		'''
		current = os.getcwd()
		pdbfilelist = os.listdir(directory)
		io = Bio.PDB.PDBIO()
		os.chdir(directory)
		print('\x1b[32m[+] Extracting files\x1b[0m')
		for TheFile in tqdm.tqdm(pdbfilelist):
			try:
				TheName = TheFile.split('.')[0].split('pdb')[1].upper()
				InFile = gzip.open(TheFile, 'rt')
				structure = Bio.PDB.PDBParser(QUIET=True).get_structure(TheName, InFile)
				count = 0
				for chain in structure.get_chains():
					io.set_structure(chain)
					io.save(structure.get_id()+'_'+chain.get_id()+'.pdb')
				os.remove(TheFile)
			except Exception as TheError:
				print('\x1b[31m[-] Failed to extract\t{}\x1b[33m{}\x1b[0m'.format(TheFile.upper(), str(TheError)))
				os.remove(TheFile)
		os.chdir(current)
	def NonProtein(self, directory):
		''' Remove non-protein structures '''
		current = os.getcwd()
		pdbfilelist = os.listdir(directory)
		os.chdir(directory)
		print('\x1b[32m[+] Deleting none-protein structures\x1b[0m')
		for TheFile in tqdm.tqdm(pdbfilelist):
			structure = Bio.PDB.PDBParser(QUIET=True).get_structure('X', TheFile)
			ppb = Bio.PDB.Polypeptide.PPBuilder()
			Type = ppb.build_peptides(structure, aa_only=True)
			if Type == []:
				os.remove(TheFile)
			else:
				continue
		os.chdir(current)
	def Break(self, directory):
		''' Remove structures with a broken (non-continuous) chains '''
		current = os.getcwd()
		pdbfilelist = os.listdir(directory)
		os.chdir(directory)
		print('\x1b[32m[+] Removing structures with non-continuous chains\x1b[0m')
		for TheFile in tqdm.tqdm(pdbfilelist):
			structure = Bio.PDB.PDBParser(QUIET=True).get_structure('X', TheFile)
			ppb = Bio.PDB.Polypeptide.PPBuilder()
			Type = ppb.build_peptides(structure, aa_only=True)
			try:
				x = Type[1]
				os.remove(TheFile)
			except:
				continue
		os.chdir(current)
	def Renumber(self, directory):
		''' Renumber structures starting at 1 '''
		current = os.getcwd()
		pdbfilelist = os.listdir(directory)
		os.chdir(directory)
		print('\x1b[32m[+] Renumbering structures\x1b[0m')
		for TheFile in tqdm.tqdm(pdbfilelist):
			pdb = open(TheFile , 'r')
			PDB = open(TheFile + 'X' , 'w')
			count = 0
			num = 0
			AA2 = None
			for line in pdb:
				count += 1
				AA1 = line[23:27]
				if not AA1 == AA2:
					num += 1
				final_line = line[:7] + '{:4d}'.format(count) + line[11:17] + line[17:21] + 'A' + '{:4d}'.format(num) + line[26:]
				AA2 = AA1
				PDB.write(final_line)
			PDB.close()
			os.remove(TheFile)
			os.rename(TheFile + 'X' , TheFile)
		os.chdir(current)

def DatasetAsPSaM(directory):
	'''
	Compile a dataset of each residue's amino acid identify, secondary
	structure, phi angle, psi angle, solvent accessible surface area as
	a .csv file and the contact map as a separate .csv file. to be run
	after clean() on the ./cleaned directory and identifying the number
	of residuis of the largest structure
	'''
	os.makedirs('./Completed', exist_ok=True)
	os.makedirs('./Error_NotEqual', exist_ok=True)
	os.makedirs('./Error_Broken', exist_ok=True)
	os.makedirs('./Error_Small', exist_ok=True)
	for File in tqdm.tqdm(os.listdir(directory)):
		try:
			TheFile = '{}/{}'.format(directory, File)
			pose = pose_from_pdb(TheFile)
			DSSP = pyrosetta.rosetta.protocols.moves.DsspMover()
			DSSP.apply(pose)
			sasa_calc = pyrosetta.rosetta.core.scoring.sasa.SasaCalc()
			sasa_calc.calculate(pose)
			size = pose.total_residue()
			aa   = []
			ss   = []
			phi  = []
			psi  = []
			sasa = []
			info = []
			ctmp = []
			m    = []
			surf = list(sasa_calc.get_residue_sasa())
			for r  in range(size):
				if pose.residue(r+1).is_protein():
					aa.append(pose.sequence(r+1, r+1))
					ss.append(pose.secstruct(r+1))
					p = pose.phi(r+1)
					if p < 0: p = p + 360
					phi.append(p)
					s = pose.psi(r+1)
					if s < 0: s = s + 360
					psi.append(s)
					sasa.append(surf[r])
			for r  in range(0, size):
				for R  in range(0, size):
					if	pose.residue(r+1).is_protein() and\
						pose.residue(R+1).is_protein():
						CAr = pose.residue(r+1).xyz('CA')
						CAR = pose.residue(R+1).xyz('CA')
						CAr_CAR_vector = CAR-CAr
						Cont = CAr_CAR_vector.norm()
						if Cont <= 12: ctmp.append(Cont)
						else: ctmp.append(0)
			if len(aa) >= 50:
					try:
						assert	len(aa) == len(ss) == len(phi)\
						== len(psi) == len(sasa) == math.sqrt(len(ctmp))
						for AA,SS,P,S,SASA in zip(aa,ss,phi,psi,sasa):
							info.append('{},{},{},{},{}'\
							.format(AA, SS, P, S, SASA))
						Info = ','.join(info)
						with open('./AsPSa.csv', 'a') as data:
							data.write(File + ',' + Info + '\n')
						with open('lengths.txt', 'a') as length:
							length.write(str(len(aa))+'\n')
						for x in ctmp:
							m.append('{}'.format(x))
						M = ','.join(m)
						with open('./M.csv', 'a') as data:
							data.write(File + ',' + M + '\n')
						os.system('mv {} ./Completed'.format(TheFile))
					except: os.system('mv {} ./Error_NotEqual'.format(TheFile))
			else: os.system('mv {} ./Error_Small'.format(TheFile))
		except: os.system('mv {} ./Error_Broken'.format(TheFile))

def Fill(filename):
    ''' Fills missing .csv table spaces with zeros '''
    with open(filename) as f:
	    with open(filename, 'a') as F:
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

def Header(length=745, choice='AsPSa'):
	'''
	Constructs a .csv header and completes the dataset. To find the value of
	the largest structure run: sort -nk 1 lengths.txt
	'''
	header = ['PDB_ID']
	if choice == 'AsPSa':
		for i in range(1, length+1):
			header.append(',aa_{},ss_{},phi_{},psi_{},sasa_{}'\
			.format(i, i, i, i, i))
		header = ''.join(header)
		with open('./AsPSa.csv', 'r') as data:
			with open('./dataset_AsPSa.csv', 'w') as head:
				head.write(header+'\n')
				for line in data:
					head.write(line)
	elif choice == 'M':
		for r in range(1, length+1):
			for c in range(1, length+1):
				header.append(',aa{}_aa{}'.format(r, c))
		header = ''.join(header)
		with open('./M.csv', 'r') as data:
			with open('./dataset_M.csv', 'w') as head:
				head.write(header+'\n')
				for line in data:
					head.write(line)
#D = Dataset()
#D.Database('DATABASE', 'PDBDatabase')
#D.Extract('PDBDatabase')
#D.NonProtein('PDBDatabase')
#D.Break('PDBDatabase')
#D.Renumber('PDBDatabase')
#DatasetAsPSaM('cln')
#Header(280, 'AsPSa')
#Fill('dataset_header.csv')
