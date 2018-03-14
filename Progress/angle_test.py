import os
from pyrosetta import *
from pyrosetta.toolbox import *
init()

pose = pose_from_pdb('3hz7.pdb')

tor = list()
for aa in range(len(pose)):
	phi = pose.phi(aa + 1)
	psi = pose.psi(aa + 1)
	tor.append(phi)
	tor.append(psi)
#print(tor)

chain = list()
for aa in range(int(len(tor) / 2)):
	chain.append('V')
chain = ''.join(chain)
#print(chain)

phi = list()
psi = list()
count = 0
for x in range(int(len(tor) / 2)):
	phi.append(tor[count])
	count += 1
	psi.append(tor[count])
	count += 1
#print(phi , psi)

pose2 = pose_from_sequence(chain)
#print(pose2)

count = 1
for aa in zip(phi , psi):
	pose2.set_phi(count , aa[0])
	pose2.set_psi(count , aa[1])
	count += 1
pose2.dump_pdb('X.pdb')

#Check
for aa in range(len(pose)):
	phi1 = pose.phi(aa + 1)
	phi2 = pose2.phi(aa + 1)

	psi1 = pose.psi(aa + 1)
	psi2 = pose2.psi(aa + 1)

	if phi1 == phi2 and psi1 == psi2:
		print('True')
#		print('res' , aa + 1 , ':\t' , round(phi1 , 3) ,'=', round(phi2 , 3) , '\t' , round(psi1 , 3) ,'=', round(psi2 , 3))

os.system('pymol X.pdb')
os.system('rm X.pdb')
