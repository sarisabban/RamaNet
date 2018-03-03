import numpy ; import matplotlib.pyplot as plt ; from mpl_toolkits.mplot3d import Axes3D
#Initial and Terminus Coordinates
O = numpy.array([1 , 1 , 2])
P = numpy.array([2 , 3 , 5])

ax = plt.figure().add_subplot(111 , projection = '3d')
ax.plot([O[0] , P[0]] , [O[1] , P[1]] , [O[2] , P[2]] , marker = 'o')

Mag = numpy.sqrt(((P[0] - O[0])**2) + ((P[1] - O[1])**2) + ((P[2] - O[2])**2))	# 3.7416573867739413
Com = [P[0] - O[0] , P[1] - O[1] , P[2] - O[2]]					# [1, 2, 3]

#1 - Move To Axis (0 , 0 , 0)
ori = P - O
ax.plot([0 , ori[0]] , [0 , ori[1]] , [0 , ori[2]] , marker = 'o')

#2 - Define Rotation Matrix
angle = numpy.radians(90)
R = [	[numpy.cos(angle)	,	-numpy.sin(angle)	,	0] , 
	[numpy.sin(angle)	, 	 numpy.cos(angle)	,	0] ,
	[	0		,		0		,	1]	]

#3 - Rotate Matrix on XY Axis, Keeping Z Unchanged
rot = numpy.dot(R , ori)
ax.plot([0 , rot[0]] , [0 , rot[1]] , [0 , rot[2]] , marker = 'o')
#4 - Move Rotated Vector Back To Original Start Point
back = rot + O
ax.plot([O[0] , back[0]] , [O[1] , back[1]] , [O[2] , back[2]] , marker = 'o')
#5 - Scale Vector To New Magnitude
Back = numpy.array([back[0] - O[0] , back[1] - O[1] , back[2] - O[2]])	#Get components of the "back" vector
Scaled = 0.5 * Back							#Multiply by scaling factor
C = numpy.add(O , Scaled)						#Add new scaled components to the initial coordinates to get new terminus coordinates
ax.plot([O[0] , C[0]] , [O[1] , C[1]] , [O[2] , C[2]] , marker = 'o')


ax.scatter(0 , 0 , 0 , marker = 's' , color = 'black' , s = 50)
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
plt.show()
