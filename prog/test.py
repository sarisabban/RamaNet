import numpy ; import matplotlib.pyplot as plt ; from mpl_toolkits.mplot3d import Axes3D

#Initial and Terminus Coordinates
O = numpy.array([1.458 , 0.000 , 0.000])
P = numpy.array([3.988 , 2.839 , 0.000])

ax = plt.figure().add_subplot(111 , projection = '3d')
ax.plot([O[0] , P[0]] , [O[1] , P[1]] , [O[2] , P[2]] , marker = 'o')

Mag = numpy.sqrt(((P[0] - O[0])**2) + ((P[1] - O[1])**2) + ((P[2] - O[2])**2))			#Magnitude = 3.8
Com = [P[0] - O[0] , P[1] - O[1] , P[2] - O[2]]

#The Carbon Atom
#1 - Move To Axis (0 , 0 , 0)
Cori = P - O
MagCori = numpy.sqrt(((Cori[0])**2) + ((Cori[1])**2) + ((Cori[2])**2))				#Magnitude = 3.8
#ax.plot([0 , Cori[0]] , [0 , Cori[1]] , [0 , Cori[2]] , marker = 'o')
#2 - Define Rotation Matrix
A = numpy.radians(00.0)	#Phi 	Angle
B = numpy.radians(20.5)	#Theta	Angle
Y = numpy.radians(00.0)	#Psi	Angle
#2D Rotation Matrix - Works Best
CR = [	[numpy.cos(B)	,	-numpy.sin(B)	,	0] , 
	[numpy.sin(B)	, 	 numpy.cos(B)	,	0] ,
	[	0	,		0	,	1]]
'''
#XYX steps Tait-Bryan angles
CR = [	[numpy.cos(B)			,	numpy.sin(B)*numpy.sin(Y)						,	numpy.sin(B)*numpy.cos(Y)					] , 
	[numpy.sin(B) *numpy.sin(A)	, 	numpy.cos(Y)*numpy.cos(A)-numpy.cos(B)*numpy.sin(Y)*numpy.sin(A)	,	numpy.cos(A)*numpy.sin(Y)-numpy.cos(B)*numpy.cos(Y)*numpy.sin(A)] ,
	[-numpy.sin(B)*numpy.cos(A)	,	numpy.cos(Y)*numpy.sin(A)+numpy.cos(B)*numpy.cos(A)*numpy.cos(Y)	,	numpy.sin(Y)*numpy.sin(A)+numpy.cos(B)*numpy.cos(Y)*numpy.cos(A)]]


#XYZ steps Tait-Bryan angles
CR = [	[numpy.cos(B)*numpy.cos(Y)						,	-numpy.cos(Y)*numpy.sin(B)						,	numpy.sin(B)			] , 
	[numpy.cos(A)*numpy.sin(Y)+numpy.cos(Y)*numpy.sin(A)*numpy.sin(B)	, 	numpy.cos(A)*numpy.cos(Y)-numpy.sin(A)*numpy.sin(B)*numpy.sin(Y)	,	-numpy.cos(B)*numpy.sin(A)	] ,
	[numpy.sin(A)*numpy.sin(Y)-numpy.cos(A)*numpy.cos(Y)*numpy.sin(B)	,	numpy.cos(Y)*numpy.sin(A)+numpy.cos(A)*numpy.sin(B)*numpy.sin(Y)	,	numpy.cos(A)*numpy.cos(B)	]]
'''
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
angle = numpy.degrees(numpy.arccos((numpy.dot(ComC , Com)) / (MagC * Mag)))			#Angle of rotation (is != 20.5 because it is calculated with the Z axis giving instead 12.08, but in 2D space the code will give 20.5)

#The Nitrogen Atom





























ax.scatter(0 , 0 , 0 , marker = 's' , color = 'black' , s = 50)
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
#plt.show()


'''
#Study!!!
#Find the A B Y angles
O =	[1.458 , 0.000 , 0.000]
C =	[2.009 , 1.420 , 0.000]
N =	[3.332 , 1.536 , 0.000]
P =	[3.988 , 2.839 , 0.000]
#The Carbon Atom
MagOP = numpy.sqrt(((P[0] - O[0])**2) + ((P[1] - O[1])**2) + ((P[2] - O[2])**2))			#Magnitude = 3.8
MagOC = numpy.sqrt(((C[0] - O[0])**2) + ((C[1] - O[1])**2) + ((C[2] - O[2])**2))			#Magnitude = 1.5
ComOP = [P[0] - O[0] , P[1] - O[1] , P[2] - O[2]]
ComOC = [C[0] - O[0] , C[1] - O[1] , C[2] - O[2]]
angle = numpy.degrees(numpy.arccos((numpy.dot(ComOP , ComOC)) / (MagOC * MagOP)))			
print(angle)
ax.plot([O[0] , P[0]] , [O[1] , P[1]] , [O[2] , P[2]] , marker = 'o')
ax.plot([O[0] , C[0]] , [O[1] , C[1]] , [O[2] , C[2]] , marker = 'o')
#The Nitrogen Atom
MagON = numpy.sqrt(((N[0] - O[0])**2) + ((N[1] - O[1])**2) + ((N[2] - O[2])**2))			#Magnitude = 2.4
ComON = [N[0] - O[0] , N[1] - O[1] , N[2] - O[2]]
angle = numpy.degrees(numpy.arccos((numpy.dot(ComOP , ComON)) / (MagON * MagOP)))			
print(angle)
ax.plot([O[0] , P[0]] , [O[1] , P[1]] , [O[2] , P[2]] , marker = 'o')
ax.plot([O[0] , N[0]] , [O[1] , N[1]] , [O[2] , N[2]] , marker = 'o')
plt.show()
'''
