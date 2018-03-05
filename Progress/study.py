import numpy , os ; import matplotlib.pyplot as plt ; from mpl_toolkits.mplot3d import Axes3D
#Initial and Terminus coordinates
O = numpy.array([3.988 , 2.839 , 0.000])
P = numpy.array([7.660 , 3.823 , 0.000])
ax = plt.figure().add_subplot(111 , projection = '3d')
ax.plot([O[0] , P[0]] , [O[1] , P[1]] , [O[2] , P[2]] , marker = 'o')
Mag = numpy.sqrt(((P[0] - O[0])**2) + ((P[1] - O[1])**2) + ((P[2] - O[2])**2))			#Magnitude = 3.8
Com = [P[0] - O[0] , P[1] - O[1] , P[2] - O[2]]
#The Carbon atom
#1 - Move To Axis (0 , 0 , 0)
Cori = P - O
MagCori = numpy.sqrt(((Cori[0])**2) + ((Cori[1])**2) + ((Cori[2])**2))				#Magnitude = 3.8
#ax.plot([0 , Cori[0]] , [0 , Cori[1]] , [0 , Cori[2]] , marker = 'o')
#2 - Define Rotation Matrix
A = numpy.radians(000.0)	#Phi 	Angle
B = numpy.radians(339.5)	#Theta	Angle
Y = numpy.radians(000.0)	#Psi	Angle
#2D Rotation Matrix - works best because the peptide bond is on a plane anyway, so there is no change in the Z axis.
"""
CR = [	[numpy.cos(B)	,	-numpy.sin(B)	,	0] , 
	[numpy.sin(B)	, 	 numpy.cos(B)	,	0] ,
	[	0	,		0	,	1]]
#XYX steps Tait-Bryan angles
CR = [	[numpy.cos(B)			,	numpy.sin(B)*numpy.sin(Y)						,	numpy.sin(B)*numpy.cos(Y)					] , 
	[numpy.sin(B) *numpy.sin(A)	, 	numpy.cos(Y)*numpy.cos(A)-numpy.cos(B)*numpy.sin(Y)*numpy.sin(A)	,	numpy.cos(A)*numpy.sin(Y)-numpy.cos(B)*numpy.cos(Y)*numpy.sin(A)] ,
	[-numpy.sin(B)*numpy.cos(A)	,	numpy.cos(Y)*numpy.sin(A)+numpy.cos(B)*numpy.cos(A)*numpy.cos(Y)	,	numpy.sin(Y)*numpy.sin(A)+numpy.cos(B)*numpy.cos(Y)*numpy.cos(A)]]
"""
#XYZ steps Tait-Bryan angles
CR = [	[numpy.cos(B)*numpy.cos(Y)						,	-numpy.cos(Y)*numpy.sin(B)						,	numpy.sin(B)			] , 
	[numpy.cos(A)*numpy.sin(Y)+numpy.cos(Y)*numpy.sin(A)*numpy.sin(B)	, 	numpy.cos(A)*numpy.cos(Y)-numpy.sin(A)*numpy.sin(B)*numpy.sin(Y)	,	-numpy.cos(B)*numpy.sin(A)	] ,
	[numpy.sin(A)*numpy.sin(Y)-numpy.cos(A)*numpy.cos(Y)*numpy.sin(B)	,	numpy.cos(Y)*numpy.sin(A)+numpy.cos(A)*numpy.sin(B)*numpy.sin(Y)	,	numpy.cos(A)*numpy.cos(B)	]]

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
AA = numpy.radians(000.0)	#Phi 	Angle
BB = numpy.radians(225.0)	#Theta	Angle
YY = numpy.radians(000.0)	#Psi	Angle
#2D Rotation Matrix - works best because the peptide bond is on a plane anyway, so there is no change in the Z axis.
"""
NR = [	[numpy.cos(BB)	,	-numpy.sin(BB)	,	0] , 
	[numpy.sin(BB)	, 	 numpy.cos(BB)	,	0] ,
	[	0	,		0	,	1]]
#XYX steps Tait-Bryan angles
NR = [	[numpy.cos(BB)			,	numpy.sin(BB)*numpy.sin(YY)						,	numpy.sin(BB)*numpy.cos(YY)						] , 
	[numpy.sin(BB) *numpy.sin(AA)	, 	numpy.cos(YY)*numpy.cos(AA)-numpy.cos(BB)*numpy.sin(YY)*numpy.sin(AA)	,	numpy.cos(AA)*numpy.sin(YY)-numpy.cos(BB)*numpy.cos(YY)*numpy.sin(AA)	] ,
	[-numpy.sin(BB)*numpy.cos(AA)	,	numpy.cos(YY)*numpy.sin(AA)+numpy.cos(BB)*numpy.cos(AA)*numpy.cos(YY)	,	numpy.sin(YY)*numpy.sin(AA)+numpy.cos(BB)*numpy.cos(YY)*numpy.cos(AA)	]]
"""
#XYZ steps Tait-Bryan angles
NR = [	[numpy.cos(BB)*numpy.cos(YY)						,	-numpy.cos(YY)*numpy.sin(BB)						,	numpy.sin(BB)			] , 
	[numpy.cos(AA)*numpy.sin(YY)+numpy.cos(YY)*numpy.sin(AA)*numpy.sin(BB)	, 	numpy.cos(AA)*numpy.cos(YY)-numpy.sin(AA)*numpy.sin(BB)*numpy.sin(YY)	,	-numpy.cos(BB)*numpy.sin(AA)	] ,
	[numpy.sin(AA)*numpy.sin(YY)-numpy.cos(AA)*numpy.cos(YY)*numpy.sin(BB)	,	numpy.cos(YY)*numpy.sin(AA)+numpy.cos(AA)*numpy.sin(BB)*numpy.sin(YY)	,	numpy.cos(AA)*numpy.cos(BB)	]]

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
NScaled = (1.5 / 3.8) * NBack										#Multiply by (distance to final N position / distance between the 2 CA atoms) the scaling factor, this gives the value each axis needs to move by to each coordinates that result in the new scaled magnitude
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
A = numpy.radians(000.0)	#Phi 	Angle
B = numpy.radians(313.3)	#Theta	Angle
Y = numpy.radians(000.0)	#Psi	Angle
#2D Rotation Matrix - works best because the peptide bond is on a plane anyway, so there is no change in the Z axis.
"""
OxR = [	[numpy.cos(B)	,	-numpy.sin(B)	,	0] , 
	[numpy.sin(B)	, 	 numpy.cos(B)	,	0] ,
	[	0	,		0	,	1]]
#XYX steps Tait-Bryan angles
OxR = [	[numpy.cos(B)			,	numpy.sin(B)*numpy.sin(Y)						,	numpy.sin(B)*numpy.cos(Y)					] , 
	[numpy.sin(B) *numpy.sin(A)	, 	numpy.cos(Y)*numpy.cos(A)-numpy.cos(B)*numpy.sin(Y)*numpy.sin(A)	,	numpy.cos(A)*numpy.sin(Y)-numpy.cos(B)*numpy.cos(Y)*numpy.sin(A)] ,
	[-numpy.sin(B)*numpy.cos(A)	,	numpy.cos(Y)*numpy.sin(A)+numpy.cos(B)*numpy.cos(A)*numpy.cos(Y)	,	numpy.sin(Y)*numpy.sin(A)+numpy.cos(B)*numpy.cos(Y)*numpy.cos(A)]]
"""
#XYZ steps Tait-Bryan angles
OxR = [	[numpy.cos(B)*numpy.cos(Y)						,	-numpy.cos(Y)*numpy.sin(B)						,	numpy.sin(B)			] , 
	[numpy.cos(A)*numpy.sin(Y)+numpy.cos(Y)*numpy.sin(A)*numpy.sin(B)	, 	numpy.cos(A)*numpy.cos(Y)-numpy.sin(A)*numpy.sin(B)*numpy.sin(Y)	,	-numpy.cos(B)*numpy.sin(A)	] ,
	[numpy.sin(A)*numpy.sin(Y)-numpy.cos(A)*numpy.cos(Y)*numpy.sin(B)	,	numpy.cos(Y)*numpy.sin(A)+numpy.cos(A)*numpy.sin(B)*numpy.sin(Y)	,	numpy.cos(A)*numpy.cos(B)	]]
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







#Study!!!
#Find the A B Y angles
O  =	[3.988 , 2.839 , 0.000]
C  =	[5.504 , 2.693 , 0.000]
Oxg=	[6.030 , 1.580 , 0.000]
N  =	[3.332 , 1.536 , 0.000]
P  =	[7.660 , 3.823 , 0.000]
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
angle = numpy.degrees(numpy.arccos((numpy.dot(ComOP , ComON)) / (MagON * MagOP)))			#Angle = 131.7
print(angle)
ax.plot([O[0] , N[0]] , [O[1] , N[1]] , [O[2] , N[2]] , marker = 'o')
#The Oxygen Atom
MagOOx = numpy.sqrt(((Oxg[0] - O[0])**2) + ((Oxg[1] - O[1])**2) + ((Oxg[2] - O[2])**2))			#Magnitude = 2.4
ComOOx = [Oxg[0] - O[0] , Oxg[1] - O[1] , Oxg[2] - O[2]]
angle = numpy.degrees(numpy.arccos((numpy.dot(ComOP , ComOOx)) / (MagOOx * MagOP)))			#Angle = 46.7
print(angle)
ax.plot([O[0] , Oxg[0]] , [O[1] , Oxg[1]] , [O[2] , Oxg[2]] , marker = 'o')
plt.show()
