import numpy , os ; import matplotlib.pyplot as plt ; from mpl_toolkits.mplot3d import Axes3D

def Draw(line):
	''' Draws a protein topology in glycine given each residue's only CA atom's XYZ coordinates '''
	''' Generates the Backbone.pdb file '''
	#Import the coordinates
	line = line.split(';')
	AAs = int((len(line)) / 3)
	count_x = 0
	count_y = 1
	count_z = 2
	ResCount = 1
	AtoCount = 1
	tag = '+'
	ax = plt.figure().add_subplot(111 , projection = '3d')
	for coordinates in range(AAs):
		if tag == '+':
			Ox = round(float(line[count_x]) , 3)
			Oy = round(float(line[count_y]) , 3)
			Oz = round(float(line[count_z]) , 3)
			if Ox == '0' and Oy == '0' and Oz == '0':
				continue
			count_x += 3
			count_y += 3
			count_z += 3
			AtoCount += 1
			try:
				Px = round(float(line[count_x]) , 3)
				Py = round(float(line[count_y]) , 3)
				Pz = round(float(line[count_z]) , 3)
			except:
				pass
			#Initial and Terminus coordinates
			O = numpy.array([Ox , Oy , Oz])
			P = numpy.array([Px , Py , Pz])
			ax.plot([O[0] , P[0]] , [O[1] , P[1]] , [O[2] , P[2]] , marker = 'o')
			Mag = numpy.sqrt(((P[0] - O[0])**2) + ((P[1] - O[1])**2) + ((P[2] - O[2])**2))			#Magnitude = 3.8
			Com = [P[0] - O[0] , P[1] - O[1] , P[2] - O[2]]
			#The Carbon atom
			#1 - Move To Axis (0 , 0 , 0)
			Cori = P - O
			MagCori = numpy.sqrt(((Cori[0])**2) + ((Cori[1])**2) + ((Cori[2])**2))				#Magnitude = 3.8
			#ax.plot([0 , Cori[0]] , [0 , Cori[1]] , [0 , Cori[2]] , marker = 'o')
			#2 - Define Rotation Matrix
			A = numpy.radians(00.0)	#Phi 	Angle
			B = numpy.radians(20.5)	#Theta	Angle
			Y = numpy.radians(00.0)	#Psi	Angle
			#2D Rotation Matrix - works best because the peptide bond is on a plane anyway, so there is no change in the Z axis.
			CR = [	[numpy.cos(B)	,	-numpy.sin(B)	,	0] , 
				[numpy.sin(B)	, 	 numpy.cos(B)	,	0] ,
				[	0	,		0	,	1]]
			"""
			#XYX steps Tait-Bryan angles
			CR = [	[numpy.cos(B)			,	numpy.sin(B)*numpy.sin(Y)						,	numpy.sin(B)*numpy.cos(Y)					] , 
				[numpy.sin(B) *numpy.sin(A)	, 	numpy.cos(Y)*numpy.cos(A)-numpy.cos(B)*numpy.sin(Y)*numpy.sin(A)	,	numpy.cos(A)*numpy.sin(Y)-numpy.cos(B)*numpy.cos(Y)*numpy.sin(A)] ,
				[-numpy.sin(B)*numpy.cos(A)	,	numpy.cos(Y)*numpy.sin(A)+numpy.cos(B)*numpy.cos(A)*numpy.cos(Y)	,	numpy.sin(Y)*numpy.sin(A)+numpy.cos(B)*numpy.cos(Y)*numpy.cos(A)]]
			#XYZ steps Tait-Bryan angles
			CR = [	[numpy.cos(B)*numpy.cos(Y)						,	-numpy.cos(Y)*numpy.sin(B)						,	numpy.sin(B)			] , 
				[numpy.cos(A)*numpy.sin(Y)+numpy.cos(Y)*numpy.sin(A)*numpy.sin(B)	, 	numpy.cos(A)*numpy.cos(Y)-numpy.sin(A)*numpy.sin(B)*numpy.sin(Y)	,	-numpy.cos(B)*numpy.sin(A)	] ,
				[numpy.sin(A)*numpy.sin(Y)-numpy.cos(A)*numpy.cos(Y)*numpy.sin(B)	,	numpy.cos(Y)*numpy.sin(A)+numpy.cos(A)*numpy.sin(B)*numpy.sin(Y)	,	numpy.cos(A)*numpy.cos(B)	]]
			"""
			#3 - Rotate Matrix on XY axis, keeping Z unchanged
			Crot = numpy.dot(CR , Cori)
			MagCrot = numpy.sqrt(((Cori[0])**2) + ((Cori[1])**2) + ((Cori[2])**2))				#Magnitude = 3.8
			#ax.plot([0 , Crot[0]] , [0 , Crot[1]] , [0 , Crot[2]] , marker = 'o')
			#4 - Move rotated vector back to original start point
			Cback = Crot + O
			MagCback = numpy.sqrt(((Cback[0] - O[0])**2) + ((Cback[1] - O[1])**2) + ((Cback[2] - O[2])**2))	#Magnitude = 3.8
			#ax.plot([O[0] , Cback[0]] , [O[1] , Cback[1]] , [O[2] , Cback[2]] , marker = 'o')
			#5 - Scale vector to new magnitude
			CBack = numpy.array([Cback[0] - O[0] , Cback[1] - O[1] , Cback[2] - O[2]])			#Get components of the vector after it was returned to its original starting point "back"
			CScaled = (1.5 / 3.8) * CBack									#Multiply by (distance to final C position / distance between the 2 CA atoms) the scaling factor, this gives the value each axis needs to move by to each coordinates that result in the new scaled magnitude
			C = numpy.add(O , CScaled)									#Add scaled components to the initial coordinates to get new terminus coordinates (the final coordintes of the C atom)
			ComC = [C[0] - O[0] , C[1] - O[1] , C[2] - O[2]]
			MagC = numpy.sqrt(((C[0] - O[0])**2) + ((C[1] - O[1])**2) + ((C[2] - O[2])**2))			#Magnitude = 1.5
			ax.plot([O[0] , C[0]] , [O[1] , C[1]] , [O[2] , C[2]] , marker = 'o')
			#6- Confirm angle between vectors
			#angleC = numpy.degrees(numpy.arccos((numpy.dot(ComC , Com)) / (MagC * Mag)))			#Angle of rotation (is != 20.5 because it is calculated with the Z axis giving instead 12.08, but in 2D space the code will give 20.5) - This is commented out because it results in a RuntimeWarning through NumPy, basically the value is very very small is gets rounded up to 0 making the division impposible, it is not a math problem rather more a NumPy problem since it has low resolution
			#The Nitrogen atom
			#1 - Move To Axis (0 , 0 , 0)
			Nori = P - O
			MagNori = numpy.sqrt(((Nori[0])**2) + ((Nori[1])**2) + ((Nori[2])**2))				#Magnitude = 3.8
			#ax.plot([0 , Nori[0]] , [0 , Nori[1]] , [0 , Nori[2]] , marker = 'o')
			#2 - Define Rotation Matrix
			AA = numpy.radians(00.0)	#Phi 	Angle
			BB = numpy.radians(351.0)	#Theta	Angle
			YY = numpy.radians(00.0)	#Psi	Angle
			#2D Rotation Matrix - works best because the peptide bond is on a plane anyway, so there is no change in the Z axis.
			NR = [	[numpy.cos(BB)	,	-numpy.sin(BB)	,	0] , 
				[numpy.sin(BB)	, 	 numpy.cos(BB)	,	0] ,
				[	0	,		0	,	1]]
			"""
			#XYX steps Tait-Bryan angles
			NR = [	[numpy.cos(BB)			,	numpy.sin(BB)*numpy.sin(YY)						,	numpy.sin(BB)*numpy.cos(YY)						] , 
				[numpy.sin(BB) *numpy.sin(AA)	, 	numpy.cos(YY)*numpy.cos(AA)-numpy.cos(BB)*numpy.sin(YY)*numpy.sin(AA)	,	numpy.cos(AA)*numpy.sin(YY)-numpy.cos(BB)*numpy.cos(YY)*numpy.sin(AA)	] ,
				[-numpy.sin(BB)*numpy.cos(AA)	,	numpy.cos(YY)*numpy.sin(AA)+numpy.cos(BB)*numpy.cos(AA)*numpy.cos(YY)	,	numpy.sin(YY)*numpy.sin(AA)+numpy.cos(BB)*numpy.cos(YY)*numpy.cos(AA)	]]
			#XYZ steps Tait-Bryan angles
			NR = [	[numpy.cos(BB)*numpy.cos(YY)						,	-numpy.cos(YY)*numpy.sin(BB)						,	numpy.sin(BB)			] , 
				[numpy.cos(AA)*numpy.sin(YY)+numpy.cos(YY)*numpy.sin(AA)*numpy.sin(BB)	, 	numpy.cos(AA)*numpy.cos(YY)-numpy.sin(AA)*numpy.sin(BB)*numpy.sin(YY)	,	-numpy.cos(BB)*numpy.sin(AA)	] ,
				[numpy.sin(AA)*numpy.sin(YY)-numpy.cos(AA)*numpy.cos(YY)*numpy.sin(BB)	,	numpy.cos(YY)*numpy.sin(AA)+numpy.cos(AA)*numpy.sin(BB)*numpy.sin(YY)	,	numpy.cos(AA)*numpy.cos(BB)	]]
			"""
			#3 - Rotate Matrix on XY axis, keeping Z unchanged
			Nrot = numpy.dot(NR , Nori)
			MagNrot = numpy.sqrt(((Nori[0])**2) + ((Nori[1])**2) + ((Nori[2])**2))					#Magnitude = 3.8
			#ax.plot([0 , Nrot[0]] , [0 , Nrot[1]] , [0 , Nrot[2]] , marker = 'o')
			#4 - Move rotated vector back to original start point
			Nback = Nrot + O
			MagNback = numpy.sqrt(((Nback[0] - O[0])**2) + ((Nback[1] - O[1])**2) + ((Nback[2] - O[2])**2))		#Magnitude = 3.8
			#ax.plot([O[0] , Nback[0]] , [O[1] , Nback[1]] , [O[2] , Nback[2]] , marker = 'o')
			#5 - Scale vector to new magnitude
			NBack = numpy.array([Nback[0] - O[0] , Nback[1] - O[1] , Nback[2] - O[2]])				#Get components of the vector after it was returned to its original starting point "back"
			NScaled = (2.4 / 3.8) * NBack										#Multiply by (distance to final N position / distance between the 2 CA atoms) the scaling factor, this gives the value each axis needs to move by to each coordinates that result in the new scaled magnitude
			N = numpy.add(O , NScaled)										#Add scaled components to the initial coordinates to get new terminus coordinates (the final coordintes of the N atom)
			ComN = [N[0] - O[0] , N[1] - O[1] , N[2] - O[2]]
			MagN = numpy.sqrt(((N[0] - O[0])**2) + ((N[1] - O[1])**2) + ((N[2] - O[2])**2))				#Magnitude = 2.4
			ax.plot([O[0] , N[0]] , [O[1] , N[1]] , [O[2] , N[2]] , marker = 'o')
			#6- Confirm angle between vectors
			#angleN = numpy.degrees(numpy.arccos((numpy.dot(ComN , Com)) / (MagN * Mag)))				#Angle of rotation (is != 20.5 because it is calculated with the Z axis giving instead 12.08, but in 2D space the code will give 20.5) - This is commented out because it results in a RuntimeWarning through NumPy, basically the value is very very small is gets rounded up to 0 making the division impposible, it is not a math problem rather more a NumPy problem since it has low resolution
			#The Oxygen atom
			#1 - Move To Axis (0 , 0 , 0)
			Oxori = P - O
			MagCori = numpy.sqrt(((Oxori[0])**2) + ((Oxori[1])**2) + ((Oxori[2])**2))				#Magnitude = 3.8
			#ax.plot([0 , Oxori[0]] , [0 , Oxori[1]] , [0 , Oxori[2]] , marker = 'o')
			#2 - Define Rotation Matrix
			A = numpy.radians(00.0)	#Phi 	Angle
			B = numpy.radians(46.7)	#Theta	Angle
			Y = numpy.radians(00.0)	#Psi	Angle
			#2D Rotation Matrix - works best because the peptide bond is on a plane anyway, so there is no change in the Z axis.
			OxR = [	[numpy.cos(B)	,	-numpy.sin(B)	,	0] , 
				[numpy.sin(B)	, 	 numpy.cos(B)	,	0] ,
				[	0	,		0	,	1]]
			"""
			#XYX steps Tait-Bryan angles
			OxR = [	[numpy.cos(B)			,	numpy.sin(B)*numpy.sin(Y)						,	numpy.sin(B)*numpy.cos(Y)					] , 
				[numpy.sin(B) *numpy.sin(A)	, 	numpy.cos(Y)*numpy.cos(A)-numpy.cos(B)*numpy.sin(Y)*numpy.sin(A)	,	numpy.cos(A)*numpy.sin(Y)-numpy.cos(B)*numpy.cos(Y)*numpy.sin(A)] ,
				[-numpy.sin(B)*numpy.cos(A)	,	numpy.cos(Y)*numpy.sin(A)+numpy.cos(B)*numpy.cos(A)*numpy.cos(Y)	,	numpy.sin(Y)*numpy.sin(A)+numpy.cos(B)*numpy.cos(Y)*numpy.cos(A)]]
			#XYZ steps Tait-Bryan angles
			OxR = [	[numpy.cos(B)*numpy.cos(Y)						,	-numpy.cos(Y)*numpy.sin(B)						,	numpy.sin(B)			] , 
				[numpy.cos(A)*numpy.sin(Y)+numpy.cos(Y)*numpy.sin(A)*numpy.sin(B)	, 	numpy.cos(A)*numpy.cos(Y)-numpy.sin(A)*numpy.sin(B)*numpy.sin(Y)	,	-numpy.cos(B)*numpy.sin(A)	] ,
				[numpy.sin(A)*numpy.sin(Y)-numpy.cos(A)*numpy.cos(Y)*numpy.sin(B)	,	numpy.cos(Y)*numpy.sin(A)+numpy.cos(A)*numpy.sin(B)*numpy.sin(Y)	,	numpy.cos(A)*numpy.cos(B)	]]
			"""
			#3 - Rotate Matrix on XY axis, keeping Z unchanged
			Oxrot = numpy.dot(OxR , Oxori)
			MagOxrot = numpy.sqrt(((Oxori[0])**2) + ((Oxori[1])**2) + ((Oxori[2])**2))				#Magnitude = 3.8
			#ax.plot([0 , Oxrot[0]] , [0 , Oxrot[1]] , [0 , Oxrot[2]] , marker = 'o')
			#4 - Move rotated vector back to original start point
			Oxback = Oxrot + O
			MagOxback = numpy.sqrt(((Oxback[0] - O[0])**2) + ((Oxback[1] - O[1])**2) + ((Oxback[2] - O[2])**2))	#Magnitude = 3.8
			#ax.plot([O[0] , Oxback[0]] , [O[1] , Oxback[1]] , [O[2] , Oxback[2]] , marker = 'o')
			#5 - Scale vector to new magnitude
			OxBack = numpy.array([Oxback[0] - O[0] , Oxback[1] - O[1] , Oxback[2] - O[2]])				#Get components of the vector after it was returned to its original starting point "back"
			OxScaled = (2.4 / 3.8) * OxBack										#Multiply by (distance to final O position / distance between the 2 CA atoms) the scaling factor, this gives the value each axis needs to move by to each coordinates that result in the new scaled magnitude
			Oxg = numpy.add(O , OxScaled)										#Add scaled components to the initial coordinates to get new terminus coordinates (the final coordintes of the O atom)
			ComOx = [Oxg[0] - O[0] , Oxg[1] - O[1] , Oxg[2] - O[2]]
			MagOx = numpy.sqrt(((Oxg[0] - O[0])**2) + ((Oxg[1] - O[1])**2) + ((Oxg[2] - O[2])**2))			#Magnitude = 2.4
			ax.plot([O[0] , Oxg[0]] , [O[1] , Oxg[1]] , [O[2] , Oxg[2]] , marker = 'o')
			#6- Confirm angle between vectors
			#angleOx = numpy.degrees(numpy.arccos((numpy.dot(ComOx , Com)) / (MagOx * Mag)))			#Angle of rotation (is != 20.5 because it is calculated with the Z axis giving instead 12.08, but in 2D space the code will give 20.5) - This is commented out because it results in a RuntimeWarning through NumPy, basically the value is very very small is gets rounded up to 0 making the division impposible, it is not a math problem rather more a NumPy problem since it has low resolution
			#Plot to visualise
			ax.scatter(0 , 0 , 0 , marker = 's' , color = 'black' , s = 50)
			ax.set_xlabel('X')
			ax.set_ylabel('Y')
			ax.set_zlabel('Z')
			#plt.show()
			'''
			#Study!!!
			#Find the A B Y angles
			O  =	[1.458 , 0.000 , 0.000]
			C  =	[2.009 , 1.420 , 0.000]
			Oxg=	[1.251 , 2.390 , 0.000]
			N  =	[3.332 , 1.536 , 0.000]
			P  =	[3.988 , 2.839 , 0.000]
			#The Carbon Atom
			MagOP = numpy.sqrt(((P[0] - O[0])**2) + ((P[1] - O[1])**2) + ((P[2] - O[2])**2))			#Magnitude = 3.8
			MagOC = numpy.sqrt(((C[0] - O[0])**2) + ((C[1] - O[1])**2) + ((C[2] - O[2])**2))			#Magnitude = 1.5
			ComOP = [P[0] - O[0] , P[1] - O[1] , P[2] - O[2]]
			ComOC = [C[0] - O[0] , C[1] - O[1] , C[2] - O[2]]
			angle = numpy.degrees(numpy.arccos((numpy.dot(ComOP , ComOC)) / (MagOC * MagOP)))			#Angle = 20.5			
			print(angle)
			ax.plot([O[0] , P[0]] , [O[1] , P[1]] , [O[2] , P[2]] , marker = 'o')
			ax.plot([O[0] , C[0]] , [O[1] , C[1]] , [O[2] , C[2]] , marker = 'o')
			#The Nitrogen Atom
			MagON = numpy.sqrt(((N[0] - O[0])**2) + ((N[1] - O[1])**2) + ((N[2] - O[2])**2))			#Magnitude = 2.4
			ComON = [N[0] - O[0] , N[1] - O[1] , N[2] - O[2]]
			angle = numpy.degrees(numpy.arccos((numpy.dot(ComOP , ComON)) / (MagON * MagOP)))			#Angle = 9.0 (360 - 9 = 351)			
			print(angle)
			ax.plot([O[0] , N[0]] , [O[1] , N[1]] , [O[2] , N[2]] , marker = 'o')
			#The Oxygen Atom
			MagOOx = numpy.sqrt(((Ox[0] - O[0])**2) + ((Ox[1] - O[1])**2) + ((Ox[2] - O[2])**2))			#Magnitude = 2.4
			ComOOx = [Ox[0] - O[0] , Ox[1] - O[1] , Ox[2] - O[2]]
			angle = numpy.degrees(numpy.arccos((numpy.dot(ComOP , ComOOx)) / (MagOOx * MagOP)))			#Angle = 46.7
			print(angle)
			ax.plot([O[0] , Ox[0]] , [O[1] , Ox[1]] , [O[2] , Ox[2]] , marker = 'o')
			plt.show()
			'''
			TheLineCA = '{:6}{:5d} {:4}{:1}{:3} {:1}{:4d}{:1}   {:8}{:8}{:8}{:6}{:6}          {:2}{:2}'.format('ATOM' , AtoCount , 'CA' , '' , 'GLY' , 'A' , ResCount , '' , str(round(Ox , 3)) , str(round(Oy , 3)) , str(round(Oz , 3)) , 1.0 , 0.0 , 'C' , '') + '\n'
			AtoCount += 1
			TheLineC = '{:6}{:5d} {:4}{:1}{:3} {:1}{:4d}{:1}   {:8}{:8}{:8}{:6}{:6}          {:2}{:2}'.format('ATOM' , AtoCount , 'C' , '' , 'GLY' , 'A' , ResCount , '' , str(round(C[0] , 3)) , str(round(C[1] , 3)) , str(round(C[2] , 3)) , 1.0 , 0.0 , 'C' , '') + '\n'
			AtoCount += 1
			TheLineOxg = '{:6}{:5d} {:4}{:1}{:3} {:1}{:4d}{:1}   {:8}{:8}{:8}{:6}{:6}          {:2}{:2}'.format('ATOM' , AtoCount , 'O' , '' , 'GLY' , 'A' , ResCount , '' , str(round(Oxg[0] , 3)) , str(round(Oxg[1] , 3)) , str(round(Oxg[2] , 3)) , 1.0 , 0.0 , 'O' , '') + '\n'
			AtoCount += 1
			ResCount += 1
			TheLineN = '{:6}{:5d} {:4}{:1}{:3} {:1}{:4d}{:1}   {:8}{:8}{:8}{:6}{:6}          {:2}{:2}'.format('ATOM' , AtoCount , 'N' , '' , 'GLY' , 'A' , ResCount , '' , str(round(N[0] , 3)) , str(round(N[1] , 3)) , str(round(N[2] , 3)) , 1.0 , 0.0 , 'N' , '') + '\n'
			output = open('Backbone.pdb' , 'a')
			output.write(TheLineCA)
			output.write(TheLineC)
			output.write(TheLineOxg)
			output.write(TheLineN)
			output.close()
			tag ='-'
		else:
			Ox = round(float(line[count_x]) , 3)
			Oy = round(float(line[count_y]) , 3)
			Oz = round(float(line[count_z]) , 3)
			if Ox == '0' and Oy == '0' and Oz == '0':
				continue
			count_x += 3
			count_y += 3
			count_z += 3
			AtoCount += 1
			try:
				Px = round(float(line[count_x]) , 3)
				Py = round(float(line[count_y]) , 3)
				Pz = round(float(line[count_z]) , 3)
			except:
				pass
			#Initial and Terminus coordinates
			O = numpy.array([Ox , Oy , Oz])
			P = numpy.array([Px , Py , Pz])
			ax.plot([O[0] , P[0]] , [O[1] , P[1]] , [O[2] , P[2]] , marker = 'o')
			Mag = numpy.sqrt(((P[0] - O[0])**2) + ((P[1] - O[1])**2) + ((P[2] - O[2])**2))			#Magnitude = 3.8
			Com = [P[0] - O[0] , P[1] - O[1] , P[2] - O[2]]
			#The Carbon atom
			#1 - Move To Axis (0 , 0 , 0)
			Cori = P - O
			MagCori = numpy.sqrt(((Cori[0])**2) + ((Cori[1])**2) + ((Cori[2])**2))				#Magnitude = 3.8
			#ax.plot([0 , Cori[0]] , [0 , Cori[1]] , [0 , Cori[2]] , marker = 'o')
			#2 - Define Rotation Matrix
			A = numpy.radians(00.0)	#Phi 	Angle
			B = numpy.radians(339.5)#Theta	Angle
			Y = numpy.radians(00.0)	#Psi	Angle
			#2D Rotation Matrix - works best because the peptide bond is on a plane anyway, so there is no change in the Z axis.
			CR = [	[numpy.cos(B)	,	-numpy.sin(B)	,	0] , 
				[numpy.sin(B)	, 	 numpy.cos(B)	,	0] ,
				[	0	,		0	,	1]]
			"""
			#XYX steps Tait-Bryan angles
			CR = [	[numpy.cos(B)			,	numpy.sin(B)*numpy.sin(Y)						,	numpy.sin(B)*numpy.cos(Y)					] , 
				[numpy.sin(B) *numpy.sin(A)	, 	numpy.cos(Y)*numpy.cos(A)-numpy.cos(B)*numpy.sin(Y)*numpy.sin(A)	,	numpy.cos(A)*numpy.sin(Y)-numpy.cos(B)*numpy.cos(Y)*numpy.sin(A)] ,
				[-numpy.sin(B)*numpy.cos(A)	,	numpy.cos(Y)*numpy.sin(A)+numpy.cos(B)*numpy.cos(A)*numpy.cos(Y)	,	numpy.sin(Y)*numpy.sin(A)+numpy.cos(B)*numpy.cos(Y)*numpy.cos(A)]]
			#XYZ steps Tait-Bryan angles
			CR = [	[numpy.cos(B)*numpy.cos(Y)						,	-numpy.cos(Y)*numpy.sin(B)						,	numpy.sin(B)			] , 
				[numpy.cos(A)*numpy.sin(Y)+numpy.cos(Y)*numpy.sin(A)*numpy.sin(B)	, 	numpy.cos(A)*numpy.cos(Y)-numpy.sin(A)*numpy.sin(B)*numpy.sin(Y)	,	-numpy.cos(B)*numpy.sin(A)	] ,
				[numpy.sin(A)*numpy.sin(Y)-numpy.cos(A)*numpy.cos(Y)*numpy.sin(B)	,	numpy.cos(Y)*numpy.sin(A)+numpy.cos(A)*numpy.sin(B)*numpy.sin(Y)	,	numpy.cos(A)*numpy.cos(B)	]]
			"""
			#3 - Rotate Matrix on XY axis, keeping Z unchanged
			Crot = numpy.dot(CR , Cori)
			MagCrot = numpy.sqrt(((Cori[0])**2) + ((Cori[1])**2) + ((Cori[2])**2))				#Magnitude = 3.8
			#ax.plot([0 , Crot[0]] , [0 , Crot[1]] , [0 , Crot[2]] , marker = 'o')
			#4 - Move rotated vector back to original start point
			Cback = Crot + O
			MagCback = numpy.sqrt(((Cback[0] - O[0])**2) + ((Cback[1] - O[1])**2) + ((Cback[2] - O[2])**2))	#Magnitude = 3.8
			#ax.plot([O[0] , Cback[0]] , [O[1] , Cback[1]] , [O[2] , Cback[2]] , marker = 'o')
			#5 - Scale vector to new magnitude
			CBack = numpy.array([Cback[0] - O[0] , Cback[1] - O[1] , Cback[2] - O[2]])			#Get components of the vector after it was returned to its original starting point "back"
			CScaled = (1.5 / 3.8) * CBack									#Multiply by (distance to final C position / distance between the 2 CA atoms) the scaling factor, this gives the value each axis needs to move by to each coordinates that result in the new scaled magnitude
			C = numpy.add(O , CScaled)									#Add scaled components to the initial coordinates to get new terminus coordinates (the final coordintes of the C atom)
			ComC = [C[0] - O[0] , C[1] - O[1] , C[2] - O[2]]
			MagC = numpy.sqrt(((C[0] - O[0])**2) + ((C[1] - O[1])**2) + ((C[2] - O[2])**2))			#Magnitude = 1.5
			ax.plot([O[0] , C[0]] , [O[1] , C[1]] , [O[2] , C[2]] , marker = 'o')
			#6- Confirm angle between vectors
			#angleC = numpy.degrees(numpy.arccos((numpy.dot(ComC , Com)) / (MagC * Mag)))			#Angle of rotation (is != 20.5 because it is calculated with the Z axis giving instead 12.08, but in 2D space the code will give 20.5) - This is commented out because it results in a RuntimeWarning through NumPy, basically the value is very very small is gets rounded up to 0 making the division impposible, it is not a math problem rather more a NumPy problem since it has low resolution
			#The Nitrogen atom
			#1 - Move To Axis (0 , 0 , 0)
			Nori = P - O
			MagNori = numpy.sqrt(((Nori[0])**2) + ((Nori[1])**2) + ((Nori[2])**2))				#Magnitude = 3.8
			#ax.plot([0 , Nori[0]] , [0 , Nori[1]] , [0 , Nori[2]] , marker = 'o')
			#2 - Define Rotation Matrix
			AA = numpy.radians(00.0)	#Phi 	Angle
			BB = numpy.radians(9.0)		#Theta	Angle
			YY = numpy.radians(00.0)	#Psi	Angle
			#2D Rotation Matrix - works best because the peptide bond is on a plane anyway, so there is no change in the Z axis.
			NR = [	[numpy.cos(BB)	,	-numpy.sin(BB)	,	0] , 
				[numpy.sin(BB)	, 	 numpy.cos(BB)	,	0] ,
				[	0	,		0	,	1]]
			"""
			#XYX steps Tait-Bryan angles
			NR = [	[numpy.cos(BB)			,	numpy.sin(BB)*numpy.sin(YY)						,	numpy.sin(BB)*numpy.cos(YY)						] , 
				[numpy.sin(BB) *numpy.sin(AA)	, 	numpy.cos(YY)*numpy.cos(AA)-numpy.cos(BB)*numpy.sin(YY)*numpy.sin(AA)	,	numpy.cos(AA)*numpy.sin(YY)-numpy.cos(BB)*numpy.cos(YY)*numpy.sin(AA)	] ,
				[-numpy.sin(BB)*numpy.cos(AA)	,	numpy.cos(YY)*numpy.sin(AA)+numpy.cos(BB)*numpy.cos(AA)*numpy.cos(YY)	,	numpy.sin(YY)*numpy.sin(AA)+numpy.cos(BB)*numpy.cos(YY)*numpy.cos(AA)	]]
			#XYZ steps Tait-Bryan angles
			NR = [	[numpy.cos(BB)*numpy.cos(YY)						,	-numpy.cos(YY)*numpy.sin(BB)						,	numpy.sin(BB)			] , 
				[numpy.cos(AA)*numpy.sin(YY)+numpy.cos(YY)*numpy.sin(AA)*numpy.sin(BB)	, 	numpy.cos(AA)*numpy.cos(YY)-numpy.sin(AA)*numpy.sin(BB)*numpy.sin(YY)	,	-numpy.cos(BB)*numpy.sin(AA)	] ,
				[numpy.sin(AA)*numpy.sin(YY)-numpy.cos(AA)*numpy.cos(YY)*numpy.sin(BB)	,	numpy.cos(YY)*numpy.sin(AA)+numpy.cos(AA)*numpy.sin(BB)*numpy.sin(YY)	,	numpy.cos(AA)*numpy.cos(BB)	]]
			"""
			#3 - Rotate Matrix on XY axis, keeping Z unchanged
			Nrot = numpy.dot(NR , Nori)
			MagNrot = numpy.sqrt(((Nori[0])**2) + ((Nori[1])**2) + ((Nori[2])**2))					#Magnitude = 3.8
			#ax.plot([0 , Nrot[0]] , [0 , Nrot[1]] , [0 , Nrot[2]] , marker = 'o')
			#4 - Move rotated vector back to original start point
			Nback = Nrot + O
			MagNback = numpy.sqrt(((Nback[0] - O[0])**2) + ((Nback[1] - O[1])**2) + ((Nback[2] - O[2])**2))		#Magnitude = 3.8
			#ax.plot([O[0] , Nback[0]] , [O[1] , Nback[1]] , [O[2] , Nback[2]] , marker = 'o')
			#5 - Scale vector to new magnitude
			NBack = numpy.array([Nback[0] - O[0] , Nback[1] - O[1] , Nback[2] - O[2]])				#Get components of the vector after it was returned to its original starting point "back"
			NScaled = (2.4 / 3.8) * NBack										#Multiply by (distance to final N position / distance between the 2 CA atoms) the scaling factor, this gives the value each axis needs to move by to each coordinates that result in the new scaled magnitude
			N = numpy.add(O , NScaled)										#Add scaled components to the initial coordinates to get new terminus coordinates (the final coordintes of the N atom)
			ComN = [N[0] - O[0] , N[1] - O[1] , N[2] - O[2]]
			MagN = numpy.sqrt(((N[0] - O[0])**2) + ((N[1] - O[1])**2) + ((N[2] - O[2])**2))				#Magnitude = 2.4
			ax.plot([O[0] , N[0]] , [O[1] , N[1]] , [O[2] , N[2]] , marker = 'o')
			#6- Confirm angle between vectors
			#angleN = numpy.degrees(numpy.arccos((numpy.dot(ComN , Com)) / (MagN * Mag)))				#Angle of rotation (is != 20.5 because it is calculated with the Z axis giving instead 12.08, but in 2D space the code will give 20.5) - This is commented out because it results in a RuntimeWarning through NumPy, basically the value is very very small is gets rounded up to 0 making the division impposible, it is not a math problem rather more a NumPy problem since it has low resolution
			#The Oxygen atom
			#1 - Move To Axis (0 , 0 , 0)
			Oxori = P - O
			MagCori = numpy.sqrt(((Oxori[0])**2) + ((Oxori[1])**2) + ((Oxori[2])**2))				#Magnitude = 3.8
			#ax.plot([0 , Oxori[0]] , [0 , Oxori[1]] , [0 , Oxori[2]] , marker = 'o')
			#2 - Define Rotation Matrix
			A = numpy.radians(00.0)	#Phi 	Angle
			B = numpy.radians(313.3)#Theta	Angle
			Y = numpy.radians(00.0)	#Psi	Angle
			#2D Rotation Matrix - works best because the peptide bond is on a plane anyway, so there is no change in the Z axis.
			OxR = [	[numpy.cos(B)	,	-numpy.sin(B)	,	0] , 
				[numpy.sin(B)	, 	 numpy.cos(B)	,	0] ,
				[	0	,		0	,	1]]
			"""
			#XYX steps Tait-Bryan angles
			OxR = [	[numpy.cos(B)			,	numpy.sin(B)*numpy.sin(Y)						,	numpy.sin(B)*numpy.cos(Y)					] , 
				[numpy.sin(B) *numpy.sin(A)	, 	numpy.cos(Y)*numpy.cos(A)-numpy.cos(B)*numpy.sin(Y)*numpy.sin(A)	,	numpy.cos(A)*numpy.sin(Y)-numpy.cos(B)*numpy.cos(Y)*numpy.sin(A)] ,
				[-numpy.sin(B)*numpy.cos(A)	,	numpy.cos(Y)*numpy.sin(A)+numpy.cos(B)*numpy.cos(A)*numpy.cos(Y)	,	numpy.sin(Y)*numpy.sin(A)+numpy.cos(B)*numpy.cos(Y)*numpy.cos(A)]]
			#XYZ steps Tait-Bryan angles
			OxR = [	[numpy.cos(B)*numpy.cos(Y)						,	-numpy.cos(Y)*numpy.sin(B)						,	numpy.sin(B)			] , 
				[numpy.cos(A)*numpy.sin(Y)+numpy.cos(Y)*numpy.sin(A)*numpy.sin(B)	, 	numpy.cos(A)*numpy.cos(Y)-numpy.sin(A)*numpy.sin(B)*numpy.sin(Y)	,	-numpy.cos(B)*numpy.sin(A)	] ,
				[numpy.sin(A)*numpy.sin(Y)-numpy.cos(A)*numpy.cos(Y)*numpy.sin(B)	,	numpy.cos(Y)*numpy.sin(A)+numpy.cos(A)*numpy.sin(B)*numpy.sin(Y)	,	numpy.cos(A)*numpy.cos(B)	]]
			"""
			#3 - Rotate Matrix on XY axis, keeping Z unchanged
			Oxrot = numpy.dot(OxR , Oxori)
			MagOxrot = numpy.sqrt(((Oxori[0])**2) + ((Oxori[1])**2) + ((Oxori[2])**2))				#Magnitude = 3.8
			#ax.plot([0 , Oxrot[0]] , [0 , Oxrot[1]] , [0 , Oxrot[2]] , marker = 'o')
			#4 - Move rotated vector back to original start point
			Oxback = Oxrot + O
			MagOxback = numpy.sqrt(((Oxback[0] - O[0])**2) + ((Oxback[1] - O[1])**2) + ((Oxback[2] - O[2])**2))	#Magnitude = 3.8
			#ax.plot([O[0] , Oxback[0]] , [O[1] , Oxback[1]] , [O[2] , Oxback[2]] , marker = 'o')
			#5 - Scale vector to new magnitude
			OxBack = numpy.array([Oxback[0] - O[0] , Oxback[1] - O[1] , Oxback[2] - O[2]])				#Get components of the vector after it was returned to its original starting point "back"
			OxScaled = (2.4 / 3.8) * OxBack										#Multiply by (distance to final O position / distance between the 2 CA atoms) the scaling factor, this gives the value each axis needs to move by to each coordinates that result in the new scaled magnitude
			Oxg = numpy.add(O , OxScaled)										#Add scaled components to the initial coordinates to get new terminus coordinates (the final coordintes of the O atom)
			ComOx = [Oxg[0] - O[0] , Oxg[1] - O[1] , Oxg[2] - O[2]]
			MagOx = numpy.sqrt(((Oxg[0] - O[0])**2) + ((Oxg[1] - O[1])**2) + ((Oxg[2] - O[2])**2))			#Magnitude = 2.4
			ax.plot([O[0] , Oxg[0]] , [O[1] , Oxg[1]] , [O[2] , Oxg[2]] , marker = 'o')
			#6- Confirm angle between vectors
			#angleOx = numpy.degrees(numpy.arccos((numpy.dot(ComOx , Com)) / (MagOx * Mag)))			#Angle of rotation (is != 20.5 because it is calculated with the Z axis giving instead 12.08, but in 2D space the code will give 20.5) - This is commented out because it results in a RuntimeWarning through NumPy, basically the value is very very small is gets rounded up to 0 making the division impposible, it is not a math problem rather more a NumPy problem since it has low resolution
			#Plot to visualise
			ax.scatter(0 , 0 , 0 , marker = 's' , color = 'black' , s = 50)
			ax.set_xlabel('X')
			ax.set_ylabel('Y')
			ax.set_zlabel('Z')
			#plt.show()
			'''
			#Study!!!
			#Find the A B Y angles
			O  =	[1.458 , 0.000 , 0.000]
			C  =	[2.009 , 1.420 , 0.000]
			Oxg=	[1.251 , 2.390 , 0.000]
			N  =	[3.332 , 1.536 , 0.000]
			P  =	[3.988 , 2.839 , 0.000]
			#The Carbon Atom
			MagOP = numpy.sqrt(((P[0] - O[0])**2) + ((P[1] - O[1])**2) + ((P[2] - O[2])**2))			#Magnitude = 3.8
			MagOC = numpy.sqrt(((C[0] - O[0])**2) + ((C[1] - O[1])**2) + ((C[2] - O[2])**2))			#Magnitude = 1.5
			ComOP = [P[0] - O[0] , P[1] - O[1] , P[2] - O[2]]
			ComOC = [C[0] - O[0] , C[1] - O[1] , C[2] - O[2]]
			angle = numpy.degrees(numpy.arccos((numpy.dot(ComOP , ComOC)) / (MagOC * MagOP)))			#Angle = 20.5			
			print(angle)
			ax.plot([O[0] , P[0]] , [O[1] , P[1]] , [O[2] , P[2]] , marker = 'o')
			ax.plot([O[0] , C[0]] , [O[1] , C[1]] , [O[2] , C[2]] , marker = 'o')
			#The Nitrogen Atom
			MagON = numpy.sqrt(((N[0] - O[0])**2) + ((N[1] - O[1])**2) + ((N[2] - O[2])**2))			#Magnitude = 2.4
			ComON = [N[0] - O[0] , N[1] - O[1] , N[2] - O[2]]
			angle = numpy.degrees(numpy.arccos((numpy.dot(ComOP , ComON)) / (MagON * MagOP)))			#Angle = 9.0 (360 - 9 = 351)			
			print(angle)
			ax.plot([O[0] , N[0]] , [O[1] , N[1]] , [O[2] , N[2]] , marker = 'o')
			#The Oxygen Atom
			MagOOx = numpy.sqrt(((Ox[0] - O[0])**2) + ((Ox[1] - O[1])**2) + ((Ox[2] - O[2])**2))			#Magnitude = 2.4
			ComOOx = [Ox[0] - O[0] , Ox[1] - O[1] , Ox[2] - O[2]]
			angle = numpy.degrees(numpy.arccos((numpy.dot(ComOP , ComOOx)) / (MagOOx * MagOP)))			#Angle = 46.7
			print(angle)
			ax.plot([O[0] , Ox[0]] , [O[1] , Ox[1]] , [O[2] , Ox[2]] , marker = 'o')
			plt.show()
			'''
			TheLineCA = '{:6}{:5d} {:4}{:1}{:3} {:1}{:4d}{:1}   {:8}{:8}{:8}{:6}{:6}          {:2}{:2}'.format('ATOM' , AtoCount , 'CA' , '' , 'GLY' , 'A' , ResCount , '' , str(round(Ox , 3)) , str(round(Oy , 3)) , str(round(Oz , 3)) , 1.0 , 0.0 , 'C' , '') + '\n'
			AtoCount += 1
			TheLineC = '{:6}{:5d} {:4}{:1}{:3} {:1}{:4d}{:1}   {:8}{:8}{:8}{:6}{:6}          {:2}{:2}'.format('ATOM' , AtoCount , 'C' , '' , 'GLY' , 'A' , ResCount , '' , str(round(C[0] , 3)) , str(round(C[1] , 3)) , str(round(C[2] , 3)) , 1.0 , 0.0 , 'C' , '') + '\n'
			AtoCount += 1
			TheLineOxg = '{:6}{:5d} {:4}{:1}{:3} {:1}{:4d}{:1}   {:8}{:8}{:8}{:6}{:6}          {:2}{:2}'.format('ATOM' , AtoCount , 'O' , '' , 'GLY' , 'A' , ResCount , '' , str(round(Oxg[0] , 3)) , str(round(Oxg[1] , 3)) , str(round(Oxg[2] , 3)) , 1.0 , 0.0 , 'O' , '') + '\n'
			AtoCount += 1
			ResCount += 1
			TheLineN = '{:6}{:5d} {:4}{:1}{:3} {:1}{:4d}{:1}   {:8}{:8}{:8}{:6}{:6}          {:2}{:2}'.format('ATOM' , AtoCount , 'N' , '' , 'GLY' , 'A' , ResCount , '' , str(round(N[0] , 3)) , str(round(N[1] , 3)) , str(round(N[2] , 3)) , 1.0 , 0.0 , 'N' , '') + '\n'
			output = open('Backbone.pdb' , 'a')
			output.write(TheLineCA)
			output.write(TheLineC)
			output.write(TheLineOxg)
			output.write(TheLineN)
			output.close()
			tag ='+'
	#Remove the last 3 lines
	os.system("sed -i '$ d' ./Backbone.pdb")
	os.system("sed -i '$ d' ./Backbone.pdb")
	os.system("sed -i '$ d' ./Backbone.pdb")
	#Add the TER as the last line
	Term = open('Backbone.pdb' , 'a')
	Term.write('TER')
	Term.close()
	#plt.show()



#line = '1.408;31.978;-3.772;0.217;30.077;-0.712;3.487;28.223;-0.125;5.525;31.454;-0.244;2.976;33.180;2.017;3.870;30.639;4.741;7.480;31.838;4.478;6.357;35.494;4.598;4.379;34.916;7.757;7.163;32.888;9.363;9.630;35.708;8.629;7.746;38.148;10.877;5.565;35.742;12.794;4.914;37.456;16.100;3.814;40.635;14.352;1.539;38.885;11.838;0.158;36.538;14.543;-0.765;39.434;16.833;-2.175;41.585;14.089;-4.383;38.793;12.653;-5.542;37.663;16.064;-6.435;41.282;16.946;-8.310;41.526;13.644;-10.263;38.455;14.687;-11.061;39.683;18.188;-12.083;43.151;17.003;-14.321;41.722;14.265;-15.878;39.167;16.676;-15.583;40.654;20.179;-17.776;38.012;21.793;-15.147;35.397;20.875;-13.181;36.844;23.791;-15.321;34.657;26.082;-13.281;31.626;25.036;-10.241;33.300;26.444;-9.167;34.747;29.798;-8.470;38.350;30.621;-8.970;39.883;27.235;-12.193;41.825;27.621;-12.073;45.465;28.451;-8.690;46.087;26.756;-7.462;47.968;23.724;-5.042;46.411;21.221;-2.033;48.054;22.861;-3.165;46.817;26.326;-3.489;43.262;24.972;-0.044;43.397;23.336;1.695;44.319;26.605;0.059;41.411;28.529;2.524;38.760;29.611;0.995;35.926;27.531;-0.798;37.500;24.611;1.951;37.571;22.035;3.440;34.278;23.169;0.173;32.397;22.908;-0.579;33.972;19.460;2.821;32.894;18.308;2.257;29.315;19.387;-1.076;29.301;17.536;0.482;30.634;14.350;3.556;28.419;14.633;1.222;25.417;14.472;-1.080;26.953;11.876;1.891;27.495;9.559;2.943;23.855;10.117;-0.561;22.688;9.150;-0.490;24.932;6.073;2.979;23.656;5.162;1.970;19.981;5.491;-1.086;20.572;3.286;-1.270;19.623;-0.391;-3.201;22.030;-2.600;-4.457;23.808;0.533;-5.916;20.509;1.883;-4.642;19.877;5.456;-3.671;16.470;6.774;-6.405;14.572;8.649;-3.938;13.995;11.493;-3.449;17.744;11.902;-7.255;18.218;12.004;-7.595;15.461;14.562;-4.985;17.107;16.803;-6.807;20.391;16.613;-10.161;18.743;17.524;-8.686;16.706;20.423;-6.648;19.447;22.082;-8.145;20.127;25.503;-7.789;23.904;24.824;-10.322;23.426;21.964;-12.899;21.350;23.812;-15.454;24.210;23.841;-15.148;25.071;20.116;-17.223;24.348;17.052;-16.313;24.770;13.438;-17.942;28.200;13.421;-15.197;29.537;15.652;-12.655;28.537;13.021;-14.715;29.875;10.161;-14.867;33.335;11.870;-11.087;33.217;12.449;-10.495;32.555;8.769;-12.857;35.314;7.683;-10.781;37.748;9.689;-7.500;36.307;8.409;-8.579;36.616;4.779;-9.899;40.166;5.363;-6.577;41.125;6.901;-4.651;39.649;3.959;-6.867;41.471;1.364;-6.305;44.798;3.109;-2.563;44.253;3.627;0.219;45.842;1.559;1.660;42.320;1.111;0.542;40.307;-2.083;-0.819;37.388;-0.063;-2.710;34.756;-2.018;-5.863;35.493;-0.006;-8.078;33.226;-2.138;-5.878;30.220;-1.242;-6.021;31.100;2.459;-9.816;31.234;2.114;-9.742;27.773;0.518;-7.377;26.469;3.196;-9.676;27.709;5.891;-12.687;26.044;4.210;-10.761;22.816;3.834;-9.597;22.942;7.485;-13.182;23.461;8.625;-14.343;20.481;6.612;-11.430;18.396;7.990;-12.349;19.440;11.529;-15.930;18.308;10.674;-14.581;14.987;9.468;-12.546;14.656;12.620;-15.636;15.082;14.768;-15.846;18.748;15.692;-19.478;19.937;15.903'
#line =  '1.408;31.978;-3.772;0.217;30.077;-0.712;3.487;28.223;-0.125;5.525;31.454;-0.244;2.976;33.180;2.017;3.870;30.639;4.741'
line =  '1.408;31.978;-3.772;0.217;30.077;-0.712;3.487;28.223;-0.125'
#line = '1.458;0.000;0.000;3.988;2.839;0.000;7.660;3.823;0.000;10.190;6.661;-0.000;13.862;7.646;0.000'

Draw(line)
#os.system('pymol Backbone.pdb')
os.system('rm Backbone.pdb')
