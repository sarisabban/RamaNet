#!/usr/bin/python3

#Isolate XYZ coordinates of CA atoms
data = open('test.pdb' , 'r')
XYZ = list()
for line in data:
	try:
		line = line.strip()
		line = line.split()
		if line[2] == 'CA':
			#Identify CA XYZ
			xyz = (float(line[6]) , float(line[7]) , float(line[8]))
			XYZ.append(xyz)
	except:
			continue

#Draw using XYZ coordinates of the CA atoms of each residue
ResCount = 1
AtoCount = 1
for CA in XYZ:
	XN = round((CA[0] - 1.411) , 3)
	YN = round((CA[1] - 0.2) , 3)
	ZN = round((CA[2] + 0.194) , 3)
	lineN = '{:6}{:5d} {:4}{:1}{:3} {:1}{:4d}{:1}   {:8}{:8}{:8}{:6}{:6}          {:2}{:2}'.format('ATOM' , AtoCount , 'N' , '' , 'GLY' , 'A' , ResCount , '' , XN , YN , ZN , 1.0 , 0.0 , 'N' , '') + '\n'
	AtoCount += 1
	lineCA = '{:6}{:5d} {:4}{:1}{:3} {:1}{:4d}{:1}   {:8}{:8}{:8}{:6}{:6}          {:2}{:2}'.format('ATOM' , AtoCount , 'CA' , '' , 'GLY' , 'A' , ResCount , '' , CA[0] , CA[1] , CA[2] , 1.0 , 0.0 , 'C' , '') + '\n'
	AtoCount += 1
	XC = round((CA[0] + 0.583) , 3)
	YC = round((CA[1] + 1.218) , 3)
	ZC = round((CA[2] + 0.716) , 3)
	lineC = '{:6}{:5d} {:4}{:1}{:3} {:1}{:4d}{:1}   {:8}{:8}{:8}{:6}{:6}          {:2}{:2}'.format('ATOM' , AtoCount , 'C' , '' , 'GLY' , 'A' , ResCount , '' , XC , YC , ZC , 1.0 , 0.0 , 'C' , '') + '\n'
	AtoCount += 1
	XO = round((CA[0] - 0.062) , 3)
	YO = round((CA[1] + 1.874) , 3)
	ZO = round((CA[2] + 1.577) , 3)
	lineO = '{:6}{:5d} {:4}{:1}{:3} {:1}{:4d}{:1}   {:8}{:8}{:8}{:6}{:6}          {:2}{:2}'.format('ATOM' , AtoCount , 'O' , '' , 'GLY' , 'A' , ResCount , '' , XO , YO , ZO , 1.0 , 0.0 , 'O' , '') + '\n'
	AtoCount += 1
	ResCount += 1
	output = open('out.pdb' , 'a')
	output.write(lineN)
	output.write(lineCA)
	output.write(lineC)
	output.write(lineO)
	output.close()

#draw glycine normally, then rotate to allow all residues to touch each other
