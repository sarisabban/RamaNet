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



	PyMol = open('PyMolMutate.py' , 'w')
	ThePyMolLine = "import pymol\ncmd.wizard('mutagenesis')\ncmd.load('Backbone.pdb')\ncmd.refresh_wizard()\nfor res in range(138): cmd.get_wizard().do_select('/Backbone//A/GLY`{}'.format(res)) ; cmd.get_wizard().set_mode('VAL') ; cmd.get_wizard().apply()\ncmd.set_wizard()\ncmd.save('Topology.pdb')"
	PyMol.write(ThePyMolLine)
	PyMol.close()
	os.system('pymol -c PyMolMutate.py')
	os.system('rm PyMolMutate.py')
	os.system('rm Backbone.pdb')
	





#Tests of 10 structures
#line = '-3.105;58.674;238.506;-1.143;56.752;235.792;-4.479;56.214;233.910;-4.885;60.080;233.914;-1.242;60.547;232.768;-1.745;58.188;229.772;-5.107;59.856;228.902;-3.380;63.295;228.842;-0.538;61.888;226.635;-3.074;60.230;224.227;-4.984;63.550;223.985;-1.770;65.511;223.120;-0.713;62.825;220.600;-4.146;63.022;218.842;-3.978;66.860;218.702;-0.384;66.706;217.268;-1.488;64.216;214.586;-4.472;66.466;213.602;-2.324;69.622;213.546;0.415;68.045;211.383;-2.262;66.639;209.027;-3.841;70.099;208.561;-0.512;71.998;208.204;0.738;69.457;205.616;-2.533;69.610;203.654;-2.233;73.437;203.354;1.295;73.042;201.837;-0.118;70.334;199.486;-2.720;72.881;198.339;0.106;75.302;197.331;1.588;72.562;195.024;-1.872;71.868;193.516;-2.498;75.589;192.967;0.845;76.095;191.211;0.224;73.048;188.993;-3.364;74.268;188.204;-2.095;77.754;187.209;0.734;76.375;185.019;-1.762;74.088;183.174;-4.357;76.888;182.837;-1.681;79.170;181.284;-0.701;76.424;178.819;-4.363;75.780;177.894;-5.000;79.523;177.314;-1.978;79.916;174.982;-2.687;76.594;173.202;-6.188;77.868;172.338;-4.671;81.130;170.916;-2.029;79.205;168.899;-4.749;76.804;167.592;-6.834;79.812;166.401;-3.735;81.213;164.538;-2.856;77.861;162.928;-6.522;77.254;161.873;-6.805;80.775;160.374;-3.478;80.228;158.542;-4.694;76.880;157.106;-7.984;78.480;155.876;-6.028;81.380;154.281;-3.825;78.857;152.430;-6.970;77.066;151.132;-8.357;80.433;149.846;-4.914;81.225;148.216;-4.952;77.785;146.465;-8.505;78.409;145.091;-7.582;81.894;143.763;-4.382;80.576;142.085;-6.550;77.817;140.441;-8.845;80.581;139.052;-5.904;82.603;137.555;-4.321;79.502;135.905;-7.748;78.761;134.260;-8.308;82.423;133.178;-4.744;82.550;131.732;-5.145;79.271;129.715;-8.269;80.739;127.996;-6.572;84.160;127.365;-3.334;82.445;126.128;-5.388;80.577;123.420;-7.064;83.865;122.286;-3.653;85.706;122.164;-2.028;82.855;120.144;-5.002;82.695;117.713;-4.879;86.503;117.044;-1.176;86.122;116.053;-2.048;83.122;113.770;-4.859;85.188;112.102;-2.379;88.087;111.456;0.096;85.692;109.781;-2.799;84.321;107.616;-3.831;87.938;106.746;-0.239;88.708;105.590;-0.296;85.492;103.456;-3.755;86.512;101.990;-2.576;90.058;101.057;0.717;88.731;99.580;-1.257;86.154;97.501;-3.842;88.742;96.367;-0.977;91.054;95.291;0.655;88.241;93.232;-2.770;87.367;91.715;-3.286;91.116;90.841;0.116;91.330;89.063;-0.606;88.030;87.229;-4.101;89.259;86.130;-2.511;92.468;84.694;0.043;90.287;82.745;-2.722;88.030;81.318;-4.695;91.174;80.288;-1.553;92.654;78.595;-1.193;89.415;76.474;-4.845;89.444;75.343;-4.614;93.214;74.610;-1.354;92.769;72.646;-3.078;90.326;70.242;-6.196;92.473;69.878;-4.123;95.605;69.232;-1.963;93.844;66.586;-5.172;92.836;64.741;-6.578;96.394;64.966;-3.208;97.901;63.855;-3.133;95.543;60.821;-6.670;96.668;59.840;-5.630;100.324;60.232;-2.537;99.669;57.993;-4.768;98.065;55.298;-7.168;101.051;55.399;-4.319;103.588;55.031;-2.742;101.596;52.186;-6.187;101.437;50.461;-6.640;105.205;50.932;-3.237;105.980;49.325;-3.888;103.480;46.486;-7.267;105.279;45.841;-5.542;108.704;45.643;-2.874;107.337;43.293;-5.267;105.504;40.910;-7.767;108.387;40.746;-5.132;111.083;40.040;-2.571;109.261;37.855;-2.274;109.654;34.053;-3.462;106.514;32.145;-1.230;104.430;29.772'
#line = '13.908;11.687;0.512;14.935;12.090;4.121;11.902;13.648;5.765;11.099;10.281;7.557;14.589;9.750;8.962;13.350;10.836;12.420;11.544;7.465;12.498;14.833;5.521;12.650;17.258;5.140;15.566;20.083;7.198;14.070;23.354;6.207;15.726;22.237;2.737;16.893;24.416;0.942;14.349;27.513;2.898;15.393;26.646;2.179;19.079;26.302;-1.582;18.447;29.570;-1.588;16.486;31.521;0.445;18.958;30.156;-0.685;22.359;28.496;-3.984;21.465;25.753;-3.534;24.059;22.386;-5.187;23.493;20.217;-3.323;26.011;18.541;-1.108;23.381;15.529;-1.630;21.045;16.442;0.240;17.900;14.079;1.027;15.053;15.593;0.163;11.689;14.675;-0.011;8.037;15.287;-3.228;6.074;17.285;-2.451;2.901;17.686;-6.032;1.618;17.309;-9.642;2.642;19.512;-12.329;1.118;19.891;-16.052;1.815;23.179;-17.186;3.271;24.596;-20.010;1.144;25.030;-22.538;3.947;24.297;-25.590;1.762;21.986;-27.082;4.403;18.707;-26.619;6.290;17.249;-24.398;7.511;17.331;-21.390;5.200;19.123;-18.516;6.884;18.807;-14.883;5.977;21.358;-12.038;6.129;19.955;-8.477;6.205;21.306;-5.092;5.158;19.638;-2.658;7.495;19.646;1.123;8.043;19.285;3.321;11.096;18.982;6.451;8.905;22.639;7.393;9.531;24.425;4.176;8.558;23.803;0.822;6.981;25.344;-2.403;8.313;25.366;-5.742;6.456;26.154;-8.207;9.235;22.693;-8.992;10.586;21.429;-12.559;10.642;17.799;-13.378;11.244;17.040;-16.072;13.889;15.391;-19.229;12.529;12.122;-18.754;14.498;11.558;-15.486;12.585;11.870;-17.043;9.095;8.085;-17.058;8.542;8.186;-13.236;8.501;10.495;-13.434;5.459;8.679;-16.314;3.787;5.229;-14.628;4.152;6.690;-11.211;3.281;5.832;-9.293;6.447;9.526;-8.556;6.736;10.742;-7.235;3.350;12.819;-4.470;1.772;11.690;-0.951;2.846;9.724;-2.078;5.935;10.471;-0.441;9.296;10.910;-2.712;12.307;11.764;-2.502;15.976;14.342;-4.962;17.332;13.162;-6.356;20.673;15.985;-8.798;21.363;19.325;-9.298;19.640;22.700;-10.861;20.311;26.150;-9.488;19.376;28.832;-12.007;18.497;32.306;-10.612;17.820;34.411;-12.686;15.362;38.060;-13.397;15.829;38.970;-10.689;13.322;36.879;-8.042;15.161;33.825;-7.866;12.888;30.473;-8.105;14.739;27.566;-10.279;13.728;24.174;-9.178;15.036;21.631;-11.956;15.432;18.063;-10.671;15.450;16.098;-12.861;17.921;12.826;-10.922;18.369;11.412;-8.074;16.253;8.204;-6.344;15.320;7.310;-4.822;11.963;6.017;-1.306;12.710;5.278;-0.332;9.068;5.100;-2.615;6.073;7.231;-2.141;2.959;4.165;-0.959;1.046;3.199;1.557;3.739;6.724;3.043;3.819;7.174;2.917;0.039;10.736;4.248;-0.245;14.233;3.995;1.113;14.515;6.492;4.008;17.430;8.968;3.707;20.918;8.297;2.393;22.565;6.222;5.109;26.230;5.591;4.572;27.720;2.093;4.837;29.392;1.694;8.251;33.076;1.194;7.981;32.896;-1.732;9.307'
#line = '-7.043;28.624;-35.873;-5.385;29.834;-33.722;-8.067;28.640;-31.319;-8.264;31.458;-28.761;-8.578;33.696;-31.857;-11.521;31.843;-33.347;-13.192;32.065;-29.942;-12.209;35.742;-29.544;-13.900;36.487;-32.827;-17.032;34.593;-31.787;-17.211;36.529;-28.518;-16.865;39.785;-30.522;-19.615;38.955;-33.005;-22.013;37.916;-30.216;-21.385;41.218;-28.427;-21.586;43.276;-31.583;-25.105;41.926;-32.211;-26.444;43.290;-28.924;-28.618;46.417;-29.383;-26.724;49.455;-28.125;-23.214;48.103;-28.780;-21.094;50.305;-31.106;-18.135;48.008;-31.471;-15.945;45.916;-29.245;-12.961;48.082;-28.219;-10.849;45.521;-26.330;-11.145;41.839;-25.384;-8.912;40.251;-22.765;-9.458;36.498;-22.566;-7.539;34.061;-20.381;-7.687;30.257;-20.613;-8.249;28.540;-17.323;-7.765;24.996;-18.670;-10.582;24.699;-21.265;-12.579;27.509;-19.627;-12.435;31.138;-20.508;-12.473;34.089;-18.183;-12.919;37.381;-20.049;-12.930;41.167;-19.616;-14.371;43.290;-22.424;-14.591;46.952;-23.191;-17.444;47.941;-25.476;-18.131;51.287;-27.176;-21.580;52.852;-26.753;-23.091;56.049;-28.176;-21.892;59.122;-26.209;-25.576;60.167;-26.124;-26.241;57.323;-23.674;-26.757;58.182;-20.010;-24.740;56.238;-17.472;-27.853;54.228;-16.639;-28.236;53.276;-20.337;-24.563;52.346;-20.725;-24.816;49.993;-17.710;-27.779;48.147;-19.209;-25.910;47.630;-22.505;-23.031;46.044;-20.514;-25.480;43.731;-18.664;-27.055;42.703;-21.963;-23.675;41.947;-23.492;-22.707;39.809;-20.510;-25.947;37.773;-20.710;-25.752;37.382;-24.482;-22.075;36.315;-24.705;-22.347;33.972;-21.711;-25.479;32.427;-23.262;-23.830;31.875;-26.661;-20.211;31.008;-26.092;-19.498;27.559;-24.650;-17.038;27.562;-21.685;-17.149;31.196;-20.741;-16.743;30.767;-17.017;-16.756;34.516;-16.324;-17.182;37.734;-18.303;-16.914;41.354;-17.177;-18.115;44.087;-19.499;-17.509;47.747;-19.188;-18.357;50.642;-21.456;-16.646;53.630;-23.095;-18.415;56.454;-24.855;-17.258;56.331;-28.572;-17.154;58.091;-31.115;-18.181;60.910;-33.272;-16.372;58.756;-35.897;-13.490;56.921;-34.308;-12.813;54.294;-35.314;-11.181;51.139;-34.553;-11.161;47.353;-34.904;-11.276;45.263;-31.698;-7.940;44.668;-29.943;-7.088;41.573;-27.885;-4.804;40.587;-24.989'
#line = '1.408;31.978;-3.772;0.217;30.077;-0.712;3.487;28.223;-0.125;5.525;31.454;-0.244;2.976;33.180;2.017;3.870;30.639;4.741;7.480;31.838;4.478;6.357;35.494;4.598;4.379;34.916;7.757;7.163;32.888;9.363;9.630;35.708;8.629;7.746;38.148;10.877;5.565;35.742;12.794;4.914;37.456;16.100;3.814;40.635;14.352;1.539;38.885;11.838;0.158;36.538;14.543;-0.765;39.434;16.833;-2.175;41.585;14.089;-4.383;38.793;12.653;-5.542;37.663;16.064;-6.435;41.282;16.946;-8.310;41.526;13.644;-10.263;38.455;14.687;-11.061;39.683;18.188;-12.083;43.151;17.003;-14.321;41.722;14.265;-15.878;39.167;16.676;-15.583;40.654;20.179;-17.776;38.012;21.793;-15.147;35.397;20.875;-13.181;36.844;23.791;-15.321;34.657;26.082;-13.281;31.626;25.036;-10.241;33.300;26.444;-9.167;34.747;29.798;-8.470;38.350;30.621;-8.970;39.883;27.235;-12.193;41.825;27.621;-12.073;45.465;28.451;-8.690;46.087;26.756;-7.462;47.968;23.724;-5.042;46.411;21.221;-2.033;48.054;22.861;-3.165;46.817;26.326;-3.489;43.262;24.972;-0.044;43.397;23.336;1.695;44.319;26.605;0.059;41.411;28.529;2.524;38.760;29.611;0.995;35.926;27.531;-0.798;37.500;24.611;1.951;37.571;22.035;3.440;34.278;23.169;0.173;32.397;22.908;-0.579;33.972;19.460;2.821;32.894;18.308;2.257;29.315;19.387;-1.076;29.301;17.536;0.482;30.634;14.350;3.556;28.419;14.633;1.222;25.417;14.472;-1.080;26.953;11.876;1.891;27.495;9.559;2.943;23.855;10.117;-0.561;22.688;9.150;-0.490;24.932;6.073;2.979;23.656;5.162;1.970;19.981;5.491;-1.086;20.572;3.286;-1.270;19.623;-0.391;-3.201;22.030;-2.600;-4.457;23.808;0.533;-5.916;20.509;1.883;-4.642;19.877;5.456;-3.671;16.470;6.774;-6.405;14.572;8.649;-3.938;13.995;11.493;-3.449;17.744;11.902;-7.255;18.218;12.004;-7.595;15.461;14.562;-4.985;17.107;16.803;-6.807;20.391;16.613;-10.161;18.743;17.524;-8.686;16.706;20.423;-6.648;19.447;22.082;-8.145;20.127;25.503;-7.789;23.904;24.824;-10.322;23.426;21.964;-12.899;21.350;23.812;-15.454;24.210;23.841;-15.148;25.071;20.116;-17.223;24.348;17.052;-16.313;24.770;13.438;-17.942;28.200;13.421;-15.197;29.537;15.652;-12.655;28.537;13.021;-14.715;29.875;10.161;-14.867;33.335;11.870;-11.087;33.217;12.449;-10.495;32.555;8.769;-12.857;35.314;7.683;-10.781;37.748;9.689;-7.500;36.307;8.409;-8.579;36.616;4.779;-9.899;40.166;5.363;-6.577;41.125;6.901;-4.651;39.649;3.959;-6.867;41.471;1.364;-6.305;44.798;3.109;-2.563;44.253;3.627;0.219;45.842;1.559;1.660;42.320;1.111;0.542;40.307;-2.083;-0.819;37.388;-0.063;-2.710;34.756;-2.018;-5.863;35.493;-0.006;-8.078;33.226;-2.138;-5.878;30.220;-1.242;-6.021;31.100;2.459;-9.816;31.234;2.114;-9.742;27.773;0.518;-7.377;26.469;3.196;-9.676;27.709;5.891;-12.687;26.044;4.210;-10.761;22.816;3.834;-9.597;22.942;7.485;-13.182;23.461;8.625;-14.343;20.481;6.612;-11.430;18.396;7.990;-12.349;19.440;11.529;-15.930;18.308;10.674;-14.581;14.987;9.468;-12.546;14.656;12.620;-15.636;15.082;14.768;-15.846;18.748;15.692;-19.478;19.937;15.903'
#line = '-8.869;37.998;62.613;-6.422;40.862;61.954;-4.934;41.710;65.333;-1.804;43.632;64.256;-2.215;47.327;63.484;0.076;49.749;61.718;0.176;53.009;63.691;0.918;55.052;60.559;-2.114;54.138;58.390;-4.115;52.152;60.939;-4.453;49.123;58.650;-4.769;45.782;60.475;-3.034;42.554;59.481;-3.557;38.860;60.157;0.170;38.389;60.951;2.713;40.245;63.006;5.618;39.559;60.587;3.451;40.934;57.757;2.792;44.048;59.847;6.560;44.450;60.236;7.305;44.261;56.520;4.255;46.460;55.692;5.659;49.215;57.945;8.719;49.346;55.616;6.344;50.742;52.919;6.170;53.850;55.092;9.739;53.941;56.340;11.579;53.028;53.186;13.699;50.535;55.105;15.420;47.968;52.796;17.397;45.591;55.009;18.279;42.925;52.409;21.132;43.654;50.019;21.848;41.474;46.974;25.410;40.018;47.397;26.029;40.129;43.708;24.897;43.650;42.765;24.159;45.363;46.049;20.582;46.284;45.304;18.588;46.627;48.576;15.123;45.179;49.255;12.487;45.502;52.003;12.042;41.715;52.295;14.253;38.682;51.892;11.886;36.999;49.342;12.008;40.126;47.122;15.834;39.869;47.201;15.577;36.137;46.392;13.366;36.533;43.317;15.620;39.337;42.100;18.861;37.339;42.613;17.509;33.978;41.379;16.424;35.235;37.982;18.615;38.307;37.445;21.822;36.537;38.487;21.714;32.823;39.264;19.690;31.789;36.208;21.979;33.479;33.574;25.174;32.503;35.511;24.099;28.858;35.367;23.031;29.173;31.692;26.456;30.698;30.888;28.172;27.972;32.901;26.423;25.244;30.938;27.562;26.808;27.645;31.092;27.244;28.969;31.386;23.617;30.087;29.965;22.261;26.844;32.679;24.195;25.089;35.497;23.132;27.518;34.368;19.556;27.205;34.471;19.813;23.426;38.054;21.202;23.497;39.155;18.444;25.875;37.550;15.771;23.782;39.414;17.134;20.731;42.685;17.071;22.593;42.030;13.573;24.022;41.345;12.340;20.481;44.627;13.844;19.211;46.668;12.466;22.155;45.165;9.064;21.656;46.617;8.804;18.195;50.038;10.144;19.152;50.179;7.274;21.697;48.286;4.566;19.819;44.818;3.636;21.121'
#line = '30.691;-26.713;1.322;31.467;-27.648;-2.290;33.169;-24.380;-3.253;36.953;-24.085;-3.628;36.922;-20.347;-4.325;35.770;-17.813;-1.780;33.747;-15.081;-3.464;31.199;-12.598;-2.173;27.730;-12.167;-3.620;28.689;-8.594;-4.511;31.780;-9.835;-6.352;29.695;-12.516;-8.118;27.279;-9.848;-9.267;30.118;-7.639;-10.477;31.153;-10.383;-12.950;28.061;-9.784;-15.066;28.718;-6.089;-15.685;29.329;-5.287;-19.352;28.148;-8.660;-20.609;24.888;-10.592;-21.155;24.313;-11.347;-17.466;24.351;-7.665;-16.529;22.363;-6.747;-19.654;19.493;-9.037;-18.607;19.543;-8.022;-14.903;20.399;-11.560;-13.823;21.096;-11.213;-10.111;22.798;-13.756;-7.880;20.832;-13.798;-4.636;22.131;-16.647;-2.491;25.181;-18.847;-2.139;24.377;-21.967;-0.150;26.937;-24.481;1.022;25.992;-26.560;-2.023;23.488;-24.468;-4.010;23.599;-21.245;-6.006;20.565;-19.107;-6.765;19.899;-16.273;-9.228;16.859;-14.419;-10.548;15.994;-12.934;-13.921;13.181;-10.660;-15.063;10.447;-11.726;-17.463;10.243;-8.557;-19.580;10.911;-7.132;-23.022;14.557;-8.313;-22.883;13.660;-11.943;-22.200;10.228;-12.503;-23.788;9.454;-14.171;-27.120;7.179;-12.415;-29.686;3.985;-13.865;-28.149;4.485;-12.599;-24.598;6.010;-15.664;-22.921;9.439;-15.915;-21.355;11.736;-17.002;-24.193;12.820;-20.608;-23.663;16.319;-19.514;-24.597;16.274;-16.960;-21.784;15.699;-19.747;-19.265;18.740;-21.456;-20.711;20.811;-18.271;-20.409;19.908;-18.247;-16.693;20.652;-21.953;-16.267;23.975;-21.806;-18.112;24.889;-18.598;-16.219;24.198;-20.092;-12.807;25.641;-23.517;-13.605;28.861;-22.216;-15.181;29.540;-19.875;-12.261;28.852;-22.576;-9.663;31.168;-24.759;-11.722;33.876;-22.101;-11.702;33.488;-21.719;-7.897;34.260;-25.429;-7.560;37.621;-25.024;-9.317;40.779;-23.214;-8.207;43.402;-20.793;-9.504;43.647;-20.546;-13.276;41.304;-23.458;-14.023;38.593;-22.436;-16.480;35.088;-23.844;-16.816;33.119;-23.834;-20.079;29.659;-25.164;-20.898;29.938;-28.322;-22.978;26.410;-29.701;-23.141;23.052;-28.338;-22.000;19.452;-29.512;-22.139;16.260;-28.006;-20.772;12.612;-28.909;-20.342;10.148;-26.054;-20.678;7.328;-27.493;-18.581;5.383;-24.286;-17.808;5.364;-20.818;-19.375;6.414;-17.807;-17.336;4.698;-14.547;-18.193;5.591;-11.830;-15.674;7.677;-10.871;-12.640;11.044;-12.529;-12.095;12.001;-16.187;-12.319;14.387;-17.839;-9.858;16.967;-20.453;-10.801;18.239;-22.888;-8.179;21.398;-24.777;-9.052;22.678;-27.992;-7.414;25.488;-30.403;-8.293;24.660;-33.821;-9.707;28.286;-34.719;-10.320;31.496;-32.915;-9.424;34.483;-34.484;-11.187;38.115;-33.668;-11.878;37.611;-32.737;-15.523;33.814;-32.363;-15.815;30.746;-31.381;-13.802;27.027;-31.975;-14.142;24.693;-29.346;-12.665;20.921;-29.363;-12.223;19.176;-25.985;-12.455;15.477;-25.608;-11.718;13.724;-22.357;-12.623;10.563;-21.303;-10.770;8.085;-18.453;-11.374;7.000;-15.946;-8.699;4.331;-18.339;-7.483;6.929;-21.055;-6.926;5.818;-23.190;-9.867;8.524;-25.038;-11.817;8.832;-23.656;-15.344;12.176;-24.955;-16.557;14.306;-27.904;-15.497;17.794;-27.901;-16.994;21.028;-29.939;-16.921;24.443;-28.438;-17.747;27.575;-30.438;-18.421;30.814;-28.473;-18.111;34.414;-29.161;-19.091;37.196;-28.030;-16.779;40.160;-26.662;-18.749'
#line = '47.310;62.592;35.224;49.572;62.259;32.026;53.463;62.196;32.092;55.546;64.162;29.554;58.921;63.028;28.498;60.498;65.807;30.503;58.467;65.449;33.582;58.848;61.841;32.971;62.546;62.903;32.452;62.179;64.690;35.678;60.470;62.482;38.428;62.849;59.879;37.192;65.936;61.743;38.161;64.963;63.225;41.572;63.656;59.784;42.255;67.133;59.016;41.010;68.659;60.226;44.349;66.499;59.012;47.226;66.809;56.352;49.943;62.933;55.236;50.355;63.894;53.724;53.590;65.017;56.868;54.884;62.426;58.719;52.806;59.037;56.909;53.992;59.534;56.457;57.711;59.180;60.263;57.301;56.026;60.629;55.247;54.043;58.078;57.236;55.459;59.389;60.654;54.609;62.683;59.327;50.991;62.003;58.026;50.904;60.141;61.185;51.533;63.473;63.172;49.934;66.392;60.981;46.779;64.643;60.518;46.561;62.220;63.549;43.238;60.399;62.480;45.682;58.490;60.397;46.851;56.582;63.542;43.303;54.848;63.272;44.937;51.862;61.566;48.055;51.699;63.800;47.815;50.984;67.471;49.351;47.733;66.937;52.494;49.954;66.596;53.370;50.863;70.182;53.789;54.536;71.389;52.326;57.064;68.755;51.477;59.552;71.476;53.868;62.231;70.488;55.797;62.976;67.321;59.222;61.945;68.650;57.336;58.936;69.896;56.082;59.018;66.278;59.457;59.317;64.282;61.393;56.732;66.323;58.849;53.925;66.372;60.301;51.147;64.241;56.941;49.376;63.254;56.618;52.974;61.532;59.842;53.288;59.564;59.194;49.714;58.745;55.835;49.621;57.015;56.829;53.166;56.290;59.545;51.742;53.905;57.135;48.908;52.874;54.705;51.534;51.685;57.341;53.388;49.685;59.588;50.743;47.843;56.153;49.076;47.312;55.807;52.399;45.520;59.138;52.607;43.385;56.952;50.745;40.995;54.434;53.626;40.361;57.608;55.065;38.922;58.628;52.016;37.225;55.300;51.472;35.456;53.710;55.134;34.931;52.642;53.156;32.124;50.039;51.446;34.178;49.626;52.079;38.096;46.328;50.570;37.385;47.139;46.874;36.803;50.429;47.405;38.608;48.132;48.294;41.639;44.769;46.421;41.353;46.027;43.912;43.844;46.915;46.274;46.600;43.605;48.348;46.180;41.377;45.216;45.933;43.561;44.609;48.990;43.156;47.969;50.777;39.446;48.084;50.323;37.947;44.790;51.332;40.676;43.812;53.922;42.435;46.970;55.457;39.104;48.762;54.650;41.229;51.977;54.336;39.942;55.325;52.969;41.474;56.417;49.727;41.325;59.968;50.999;44.169;58.698;52.856;46.548;58.149;49.644;45.749;61.656;48.781;46.670;62.402;52.330;50.170;60.367;52.127;50.270;61.682;48.630;50.771;65.234;49.675;53.179;65.150;52.993;55.133;63.020;50.636;55.529;65.935;48.284;56.310;68.004;51.314;59.136;65.636;52.525;60.652;65.527;48.982;60.945;69.119;49.201;62.194;68.984;52.859;65.157;67.272;51.048;66.126;68.845;47.719;65.423;72.495;48.706;65.808;74.405;45.479;66.862;71.150;43.654;62.967;70.465;43.651;61.388;73.323;41.589;57.576;73.837;40.559;57.918;72.023;37.318;59.373;69.142;39.362;56.395;69.002;42.032;53.780;69.353;39.443;55.269;66.290;37.511;55.324;65.029;41.045;51.621;66.102;41.844;51.174;64.385;38.620;52.127;60.815;39.956;51.313;61.355;43.475;47.638;61.654;42.201;47.477;58.732;39.591;49.166;56.457;42.094;46.574;57.863;44.527;43.863;57.297;41.936;45.416;53.823;41.130;45.345;52.831;44.724;41.397;53.581;44.745;40.495;52.280;41.346;38.855;49.089;42.897;36.693;50.557;45.768;33.296;48.719;45.777'
#line = '-10.236;16.678;38.099;-7.015;15.822;36.082;-4.896;17.084;33.268;-2.529;14.895;31.169;0.220;16.138;28.826;-0.671;16.685;25.206;0.205;14.107;22.487;3.156;15.184;20.383;5.789;13.711;18.221;9.370;13.920;19.166;11.364;17.040;18.568;14.865;18.473;19.385;15.002;21.808;21.150;18.715;21.529;20.865;20.738;19.395;18.482;24.181;17.961;18.340;26.847;20.004;16.706;29.891;19.076;14.795;32.944;18.789;16.923;36.683;18.331;16.432;38.744;16.659;19.097;42.375;15.765;19.878;43.268;12.095;20.096;43.311;10.991;23.659;40.984;13.617;24.830;37.270;13.457;25.569;34.220;14.527;23.598;30.615;15.033;24.409;27.562;15.305;22.139;24.590;17.096;23.651;20.960;17.271;22.591;17.707;17.807;24.430;14.534;16.024;23.490;10.973;16.996;23.718;7.861;14.920;23.065;5.164;12.768;24.690;4.939;9.990;25.364;8.485;10.528;26.080;10.594;9.228;23.154;13.132;6.547;23.133;16.544;7.311;22.350;19.226;5.591;20.433;22.675;6.613;19.280;24.498;5.456;16.228;28.007;5.653;14.940;28.045;5.435;11.223;24.830;3.578;10.616;25.779;1.327;13.634;24.096;1.129;17.039;25.960;2.397;19.987;26.011;0.130;23.022;26.876;1.771;26.250;30.320;0.900;27.437;33.301;2.111;29.065;34.116;4.118;25.862;30.530;5.138;25.301;28.730;6.970;27.983;25.142;7.801;27.446;22.739;9.622;29.758;19.146;10.670;29.403;17.292;12.926;31.883;13.782;13.222;32.870;14.207;16.671;31.344;15.065;15.189;28.039;18.631;16.146;28.283;20.303;13.399;26.220;24.073;12.860;25.840;27.052;10.869;24.503;30.429;11.048;26.111;33.482;9.343;24.429;36.226;8.880;27.125;39.560;9.018;25.337;38.825;9.595;21.688;40.598;7.263;19.260;40.955;7.050;15.624;38.220;4.542;15.500;35.868;7.240;16.945;36.031;9.432;13.808;32.571;9.127;12.445;29.108;10.516;11.924;26.993;10.295;15.022;23.133;10.172;15.213;20.407;10.485;17.807;17.060;9.090;16.970;14.177;9.875;19.263;11.040;7.937;18.643;7.587;7.541;20.222;4.158;6.379;19.029;3.741;9.294;16.780;7.042;9.738;15.133;10.723;10.182;14.893;13.359;12.931;15.128;16.929;12.723;14.324;19.942;14.816;14.407;23.422;14.012;13.134;26.962;15.226;13.704;30.546;14.447;12.920;33.570;13.773;15.019;36.840;14.631;13.442;39.639;13.005;15.439;42.498;15.416;15.246;46.034;14.284;14.694;48.943;16.510;15.733'
#line = '97.281;48.971;20.617;98.218;45.315;21.092;99.305;45.573;24.732;96.173;43.635;25.837;97.384;40.511;24.007;99.500;37.793;25.652;101.448;37.300;22.399;102.513;40.942;22.044;106.056;41.767;23.250;106.239;45.420;22.169;109.707;46.652;21.236;111.591;43.565;22.478;114.925;43.036;20.649;116.060;39.884;18.856;118.635;38.906;21.564;116.025;38.942;24.323;113.414;37.306;22.215;115.809;34.402;21.401;116.710;33.952;25.054;113.124;34.173;26.070;111.637;31.738;23.545;114.514;29.161;23.600;113.429;28.127;27.150;109.892;27.276;26.049;108.017;24.548;24.260;107.875;24.257;20.491;105.075;26.577;19.171;105.665;29.327;21.757;105.556;32.656;19.927;105.452;36.423;20.147;105.012;39.420;17.948;107.475;42.284;18.378;107.817;45.722;16.905;110.759;46.837;14.778;110.280;50.612;15.051;113.212;51.710;12.840;112.352;49.304;9.981;110.141;49.045;6.890;108.670;45.779;8.128;107.272;47.161;11.407;106.745;43.764;13.099;108.231;40.313;13.325;106.871;37.029;14.498;109.303;34.969;16.526;108.570;31.307;17.239;110.223;28.320;18.916;110.266;25.250;16.649;112.188;22.092;17.205;115.221;23.402;19.026;115.587;26.830;17.313;113.948;30.192;17.452;112.818;31.372;14.042;111.981;34.917;12.991;109.580;36.112;10.356;109.214;39.650;8.944;105.702;41.022;8.796;105.339;43.068;5.649;102.832;45.809;4.940;102.321;46.738;1.320;99.468;48.604;-0.321;96.324;47.911;1.731;97.478;44.575;3.269;99.753;42.804;5.722;101.449;39.461;5.217;104.252;37.067;5.996;104.073;35.644;2.440;102.423;37.441;-0.496;100.193;34.410;-1.232;98.235;35.098;1.934;97.303;38.733;2.644;95.254;40.148;5.530;93.521;43.454;6.354;95.098;43.716;9.861;97.826;42.378;12.052;95.206;40.807;14.368;93.645;38.837;11.479;97.052;37.626;10.374;97.928;36.436;13.887;94.502;34.811;14.361;94.866;32.944;11.043;98.206;31.361;12.006;97.048;30.620;15.544;94.989;27.870;13.855;96.700;27.279;10.405;100.233;25.919;10.386;102.996;28.205;9.168;103.994;25.296;6.890;101.610;26.891;4.310;104.241;29.690;4.024;107.459;27.996;5.168;107.530;24.256;4.612;110.160;23.462;7.403;108.146;25.410;9.996;105.295;23.026;10.628;103.641;24.662;13.614;100.793;27.054;14.425;101.287;30.084;16.621;100.348;27.777;19.443;101.322;29.310;22.818;101.530;33.077;23.399;104.775;33.739;25.214;103.114;35.462;28.234;100.814;32.464;28.706'
#line = '28.895;30.940;32.268;26.571;33.913;31.935;27.958;36.738;29.838;26.711;40.107;28.624;27.705;41.552;25.238;26.897;45.127;24.265;27.013;45.464;20.469;26.996;48.269;17.923;26.221;46.871;14.473;26.911;48.756;11.236;26.129;48.442;7.541;22.819;48.015;5.750;21.141;47.121;9.040;17.936;49.216;9.145;14.837;48.660;7.021;15.721;44.988;6.521;13.889;43.240;9.379;17.054;43.142;11.471;15.352;43.833;14.803;12.513;41.387;14.180;14.944;38.684;13.036;17.386;39.288;15.903;14.468;39.003;18.350;13.338;35.754;16.718;16.824;34.246;16.936;17.332;35.469;20.516;14.097;33.835;21.620;15.001;30.602;19.774;18.398;30.533;21.530;17.225;31.717;24.973;19.221;34.940;24.895;17.937;37.878;26.960;18.308;41.607;26.378;17.269;44.094;23.729;17.876;45.656;20.328;17.183;49.113;18.894;17.829;51.335;15.900;20.044;54.382;16.432;19.663;57.665;14.524;23.117;57.306;12.912;21.918;54.218;10.968;23.457;51.624;13.289;21.821;49.010;15.510;22.544;48.692;19.234;21.970;45.485;21.207;22.491;43.912;24.619;22.759;40.126;24.642;22.861;38.329;27.984;22.722;34.641;28.810;24.681;31.420;28.538;27.889;32.186;26.650;27.325;29.528;23.968;23.865;30.987;23.249;25.128;34.576;23.235;27.986;33.569;20.908;25.538;31.847;18.538;23.465;35.047;18.315;26.543;37.256;17.792;27.618;35.034;14.913;24.116;35.177;13.450;24.256;38.975;13.680;27.063;39.133;11.095;24.777;37.389;8.592;21.464;38.964;9.606;18.274;38.650;7.600;18.608;38.245;3.813;18.884;41.929;2.876;21.185;42.876;5.744;24.835;43.655;4.991;26.799;43.894;8.225;30.136;45.645;7.773;31.278;45.567;11.397;30.224;44.713;14.931;31.796;46.569;17.873;31.607;44.625;21.112;32.000;45.461;24.796;31.740;42.706;27.382;29.901;43.801;30.532;28.839;42.504;33.928;25.512;40.630;33.609;22.554;42.610;35.019;19.734;40.034;35.053;18.021;40.625;31.718;14.910;38.496;31.200;13.387;37.982;27.772;14.411;39.757;24.585;12.628;42.844;23.277;12.850;45.637;20.720;13.184;49.033;22.375;13.792;51.538;19.575'

line = '1.408;31.978;-3.772;0.217;30.077;-0.712;3.487;28.223;-0.125;5.525;31.454;-0.244;2.976;33.180;2.017;3.870;30.639;4.741;7.480;31.838;4.478;6.357;35.494;4.598;4.379;34.916;7.757;7.163;32.888;9.363;9.630;35.708;8.629;7.746;38.148;10.877;5.565;35.742;12.794;4.914;37.456;16.100;3.814;40.635;14.352;1.539;38.885;11.838;0.158;36.538;14.543;-0.765;39.434;16.833;-2.175;41.585;14.089;-4.383;38.793;12.653;-5.542;37.663;16.064;-6.435;41.282;16.946;-8.310;41.526;13.644;-10.263;38.455;14.687;-11.061;39.683;18.188;-12.083;43.151;17.003;-14.321;41.722;14.265;-15.878;39.167;16.676;-15.583;40.654;20.179;-17.776;38.012;21.793;-15.147;35.397;20.875;-13.181;36.844;23.791;-15.321;34.657;26.082;-13.281;31.626;25.036;-10.241;33.300;26.444;-9.167;34.747;29.798;-8.470;38.350;30.621;-8.970;39.883;27.235;-12.193;41.825;27.621;-12.073;45.465;28.451;-8.690;46.087;26.756;-7.462;47.968;23.724;-5.042;46.411;21.221;-2.033;48.054;22.861;-3.165;46.817;26.326;-3.489;43.262;24.972;-0.044;43.397;23.336;1.695;44.319;26.605;0.059;41.411;28.529;2.524;38.760;29.611;0.995;35.926;27.531;-0.798;37.500;24.611;1.951;37.571;22.035;3.440;34.278;23.169;0.173;32.397;22.908;-0.579;33.972;19.460;2.821;32.894;18.308;2.257;29.315;19.387;-1.076;29.301;17.536;0.482;30.634;14.350;3.556;28.419;14.633;1.222;25.417;14.472;-1.080;26.953;11.876;1.891;27.495;9.559;2.943;23.855;10.117;-0.561;22.688;9.150;-0.490;24.932;6.073;2.979;23.656;5.162;1.970;19.981;5.491;-1.086;20.572;3.286;-1.270;19.623;-0.391;-3.201;22.030;-2.600;-4.457;23.808;0.533;-5.916;20.509;1.883;-4.642;19.877;5.456;-3.671;16.470;6.774;-6.405;14.572;8.649;-3.938;13.995;11.493;-3.449;17.744;11.902;-7.255;18.218;12.004;-7.595;15.461;14.562;-4.985;17.107;16.803;-6.807;20.391;16.613;-10.161;18.743;17.524;-8.686;16.706;20.423;-6.648;19.447;22.082;-8.145;20.127;25.503;-7.789;23.904;24.824;-10.322;23.426;21.964;-12.899;21.350;23.812;-15.454;24.210;23.841;-15.148;25.071;20.116;-17.223;24.348;17.052;-16.313;24.770;13.438;-17.942;28.200;13.421;-15.197;29.537;15.652;-12.655;28.537;13.021;-14.715;29.875;10.161;-14.867;33.335;11.870;-11.087;33.217;12.449;-10.495;32.555;8.769;-12.857;35.314;7.683;-10.781;37.748;9.689;-7.500;36.307;8.409;-8.579;36.616;4.779;-9.899;40.166;5.363;-6.577;41.125;6.901;-4.651;39.649;3.959;-6.867;41.471;1.364;-6.305;44.798;3.109;-2.563;44.253;3.627;0.219;45.842;1.559;1.660;42.320;1.111;0.542;40.307;-2.083;-0.819;37.388;-0.063;-2.710;34.756;-2.018;-5.863;35.493;-0.006;-8.078;33.226;-2.138;-5.878;30.220;-1.242;-6.021;31.100;2.459;-9.816;31.234;2.114;-9.742;27.773;0.518;-7.377;26.469;3.196;-9.676;27.709;5.891;-12.687;26.044;4.210;-10.761;22.816;3.834;-9.597;22.942;7.485;-13.182;23.461;8.625;-14.343;20.481;6.612;-11.430;18.396;7.990;-12.349;19.440;11.529;-15.930;18.308;10.674;-14.581;14.987;9.468;-12.546;14.656;12.620;-15.636;15.082;14.768;-15.846;18.748;15.692;-19.478;19.937;15.903'
#line =  '1.408;31.978;-3.772;0.217;30.077;-0.712;3.487;28.223;-0.125;5.525;31.454;-0.244;2.976;33.180;2.017;3.870;30.639;4.741'
#line =  '1.408;31.978;-3.772;0.217;30.077;-0.712;3.487;28.223;-0.125'
#line = '1.458;0.000;0.000;3.988;2.839;0.000;7.660;3.823;0.000;10.190;6.661;-0.000;13.862;7.646;0.000'

Draw(line)
os.system('pymol Topology.pdb')
os.system('rm Topology.pdb')
