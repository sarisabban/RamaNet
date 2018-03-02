import os , numpy
import matplotlib.pyplot as plt ; from mpl_toolkits.mplot3d import Axes3D

def DrawPDB(line):
	''' Draws a protein topology given the CA atom's XYZ coordinates of each residue '''
	''' Generates the DeNovo1.pdb file '''
	line = line.split(';')
	items = int((len(line) - 2) / 3)
	count_x = 0
	count_y = 1
	count_z = 2
	ResCount = 1
	AtoCount = 1
	for coordinates in range(items):
		x = str(round(float(line[count_x]) , 3))
		y = str(round(float(line[count_y]) , 3))
		z = str(round(float(line[count_z]) , 3))
		if x == '0' and y == '0' and z == '0':
			continue
		TheLine = '{:6}{:5d} {:4}{:1}{:3} {:1}{:4d}{:1}   {:8}{:8}{:8}{:6}{:6}          {:2}{:2}'.format('ATOM' , AtoCount , 'CA' , '' , 'GLY' , 'A' , ResCount , '' , x , y , z , 1.0 , 0.0 , 'C' , '') + '\n'
		output = open('DeNovo1.pdb' , 'a')
		output.write(TheLine)
		output.close()
		count_x += 3
		count_y += 3
		count_z += 3
		AtoCount += 1
		ResCount += 1

def ConstructPDB(filename):
	''' Uses the CA atom's XYZ coordinates of each residue to construct the atoms of a Glycine amino acid '''
	''' Generates the DeNovo2.pdb file '''
	data = open(filename , 'r')
	ResCount = 1
	AtoCount = 1
	tag = '+'
	for line in data:
		if tag == '+':
			line = line.split()
			CAx = line[6]
			CAy = line[7]
			CAz = line[8]
			Nx = round(float(line[6]) - 1.450 , 3)
			Ny = round(float(line[7]) - 0.171 , 3)
			Nz = line[8]
			Cx = round(float(line[6]) + 0.383 , 3)
			Cy = round(float(line[7]) - 1.401 , 3)
			Cz = line[8]
			Ox = round(float(line[6]) - 0.475 , 3)
			Oy = round(float(line[7]) - 2.350 , 3)
			Oz = line[8]
			TheNLine  = '{:6}{:5d} {:4}{:1}{:3} {:1}{:4d}{:1}   {:8}{:8}{:8}{:6}{:6}          {:2}{:2}'.format('ATOM' , AtoCount , 'N'  , '' , 'GLY' , 'A' , ResCount , '' , Nx  , Ny  , Nz  , 1.0 , 0.0 , 'N' , '') + '\n'
			AtoCount += 1
			TheCALine = '{:6}{:5d} {:4}{:1}{:3} {:1}{:4d}{:1}   {:8}{:8}{:8}{:6}{:6}          {:2}{:2}'.format('ATOM' , AtoCount , 'CA' , '' , 'GLY' , 'A' , ResCount , '' , CAx , CAy , CAz , 1.0 , 0.0 , 'C' , '') + '\n'
			AtoCount += 1
			TheCLine  = '{:6}{:5d} {:4}{:1}{:3} {:1}{:4d}{:1}   {:8}{:8}{:8}{:6}{:6}          {:2}{:2}'.format('ATOM' , AtoCount , 'C'  , '' , 'GLY' , 'A' , ResCount , '' , Cx  , Cy  , Cz  , 1.0 , 0.0 , 'C' , '') + '\n'
			AtoCount += 1
			TheOLine  = '{:6}{:5d} {:4}{:1}{:3} {:1}{:4d}{:1}   {:8}{:8}{:8}{:6}{:6}          {:2}{:2}'.format('ATOM' , AtoCount , 'O'  , '' , 'GLY' , 'A' , ResCount , '' , Ox  , Oy  , Oz  , 1.0 , 0.0 , 'O' , '') + '\n'
			AtoCount += 1
			ResCount += 1
			output = open('DeNovo2.pdb' , 'a')
			output.write(TheNLine)
			output.write(TheCALine)
			output.write(TheCLine)
			output.write(TheOLine)
			output.close()
			tag ='-'
		else:
			line = line.split()
			CAx = line[6]
			CAy = line[7]
			CAz = line[8]
			Nx = round(float(line[6]) + 1.450 , 3)
			Ny = round(float(line[7]) + 0.171 , 3)
			Nz = line[8]
			Cx = round(float(line[6]) - 0.383 , 3)
			Cy = round(float(line[7]) + 1.401 , 3)
			Cz = line[8]
			Ox = round(float(line[6]) + 0.475 , 3)
			Oy = round(float(line[7]) + 2.350 , 3)
			Oz = line[8]
			TheNLine  = '{:6}{:5d} {:4}{:1}{:3} {:1}{:4d}{:1}   {:8}{:8}{:8}{:6}{:6}          {:2}{:2}'.format('ATOM' , AtoCount , 'N'  , '' , 'GLY' , 'A' , ResCount , '' , Nx  , Ny  , Nz  , 1.0 , 0.0 , 'N' , '') + '\n'
			AtoCount += 1
			TheCALine = '{:6}{:5d} {:4}{:1}{:3} {:1}{:4d}{:1}   {:8}{:8}{:8}{:6}{:6}          {:2}{:2}'.format('ATOM' , AtoCount , 'CA' , '' , 'GLY' , 'A' , ResCount , '' , CAx , CAy , CAz , 1.0 , 0.0 , 'C' , '') + '\n'
			AtoCount += 1
			TheCLine  = '{:6}{:5d} {:4}{:1}{:3} {:1}{:4d}{:1}   {:8}{:8}{:8}{:6}{:6}          {:2}{:2}'.format('ATOM' , AtoCount , 'C'  , '' , 'GLY' , 'A' , ResCount , '' , Cx  , Cy  , Cz  , 1.0 , 0.0 , 'C' , '') + '\n'
			AtoCount += 1
			TheOLine  = '{:6}{:5d} {:4}{:1}{:3} {:1}{:4d}{:1}   {:8}{:8}{:8}{:6}{:6}          {:2}{:2}'.format('ATOM' , AtoCount , 'O'  , '' , 'GLY' , 'A' , ResCount , '' , Ox  , Oy  , Oz  , 1.0 , 0.0 , 'O' , '') + '\n'
			AtoCount += 1
			ResCount += 1
			output = open('DeNovo2.pdb' , 'a')
			output.write(TheNLine)
			output.write(TheCALine)
			output.write(TheCLine)
			output.write(TheOLine)
			output.close()
			tag = '+'

def RotatePDB(filename):
	''' Rotates each Glycine amino acid until they touch and form a chain '''
	''' Generates the DeNovo3.pdb file '''
	#Initial and Terminus Coordinates
	P = numpy.array([2 , 1 , 0])
	O = numpy.array([3 , 5 , 0])
	#1 - Move Point To Axis Origins
	t = O - P
	#Define Angle and Use Radians Instead of Degrees 
	angle = numpy.radians(90)
	#Define Rotation Matrix
	r = [	[numpy.cos(angle)	,	-numpy.sin(angle)	,	0] , 
		[numpy.sin(angle)	, 	 numpy.cos(angle)	,	0] ,
		[	0		,		0		,	1]	]
	#2 - Dot Product To Rotate Vector
	coo = numpy.dot(r , t)
	#3 - More Rotated Vector Back To Original Start Point
	coo2 = coo + P
	#4 - Scale Down To New Magnitude
	V = numpy.array([coo2[0] - P[0] , coo2[1] - P[1] , coo2[2] - P[2]])	#Get components of coo2
	V2 = 2 * V								#Multiply by scaling factor
	coo3 = numpy.add(P , V2)						#Add new scaled components to the initial coordinates to get new terminus coordinates




	plt.plot([P[0] , O[0]] , [P[1] , O[1]])
	plt.plot([P[0] , -6] , [P[1] , 3])
	plt.plot([P[0] , coo2[0]] , [P[1] , coo2[1]])
	plt.show()


#line = '1.408;31.978;-3.772;0.217;30.077;-0.712;3.487;28.223;-0.125;5.525;31.454;-0.244;2.976;33.180;2.017;3.870;30.639;4.741;7.480;31.838;4.478;6.357;35.494;4.598;4.379;34.916;7.757;7.163;32.888;9.363;9.630;35.708;8.629;7.746;38.148;10.877;5.565;35.742;12.794;4.914;37.456;16.100;3.814;40.635;14.352;1.539;38.885;11.838;0.158;36.538;14.543;-0.765;39.434;16.833;-2.175;41.585;14.089;-4.383;38.793;12.653;-5.542;37.663;16.064;-6.435;41.282;16.946;-8.310;41.526;13.644;-10.263;38.455;14.687;-11.061;39.683;18.188;-12.083;43.151;17.003;-14.321;41.722;14.265;-15.878;39.167;16.676;-15.583;40.654;20.179;-17.776;38.012;21.793;-15.147;35.397;20.875;-13.181;36.844;23.791;-15.321;34.657;26.082;-13.281;31.626;25.036;-10.241;33.300;26.444;-9.167;34.747;29.798;-8.470;38.350;30.621;-8.970;39.883;27.235;-12.193;41.825;27.621;-12.073;45.465;28.451;-8.690;46.087;26.756;-7.462;47.968;23.724;-5.042;46.411;21.221;-2.033;48.054;22.861;-3.165;46.817;26.326;-3.489;43.262;24.972;-0.044;43.397;23.336;1.695;44.319;26.605;0.059;41.411;28.529;2.524;38.760;29.611;0.995;35.926;27.531;-0.798;37.500;24.611;1.951;37.571;22.035;3.440;34.278;23.169;0.173;32.397;22.908;-0.579;33.972;19.460;2.821;32.894;18.308;2.257;29.315;19.387;-1.076;29.301;17.536;0.482;30.634;14.350;3.556;28.419;14.633;1.222;25.417;14.472;-1.080;26.953;11.876;1.891;27.495;9.559;2.943;23.855;10.117;-0.561;22.688;9.150;-0.490;24.932;6.073;2.979;23.656;5.162;1.970;19.981;5.491;-1.086;20.572;3.286;-1.270;19.623;-0.391;-3.201;22.030;-2.600;-4.457;23.808;0.533;-5.916;20.509;1.883;-4.642;19.877;5.456;-3.671;16.470;6.774;-6.405;14.572;8.649;-3.938;13.995;11.493;-3.449;17.744;11.902;-7.255;18.218;12.004;-7.595;15.461;14.562;-4.985;17.107;16.803;-6.807;20.391;16.613;-10.161;18.743;17.524;-8.686;16.706;20.423;-6.648;19.447;22.082;-8.145;20.127;25.503;-7.789;23.904;24.824;-10.322;23.426;21.964;-12.899;21.350;23.812;-15.454;24.210;23.841;-15.148;25.071;20.116;-17.223;24.348;17.052;-16.313;24.770;13.438;-17.942;28.200;13.421;-15.197;29.537;15.652;-12.655;28.537;13.021;-14.715;29.875;10.161;-14.867;33.335;11.870;-11.087;33.217;12.449;-10.495;32.555;8.769;-12.857;35.314;7.683;-10.781;37.748;9.689;-7.500;36.307;8.409;-8.579;36.616;4.779;-9.899;40.166;5.363;-6.577;41.125;6.901;-4.651;39.649;3.959;-6.867;41.471;1.364;-6.305;44.798;3.109;-2.563;44.253;3.627;0.219;45.842;1.559;1.660;42.320;1.111;0.542;40.307;-2.083;-0.819;37.388;-0.063;-2.710;34.756;-2.018;-5.863;35.493;-0.006;-8.078;33.226;-2.138;-5.878;30.220;-1.242;-6.021;31.100;2.459;-9.816;31.234;2.114;-9.742;27.773;0.518;-7.377;26.469;3.196;-9.676;27.709;5.891;-12.687;26.044;4.210;-10.761;22.816;3.834;-9.597;22.942;7.485;-13.182;23.461;8.625;-14.343;20.481;6.612;-11.430;18.396;7.990;-12.349;19.440;11.529;-15.930;18.308;10.674;-14.581;14.987;9.468;-12.546;14.656;12.620;-15.636;15.082;14.768;-15.846;18.748;15.692;-19.478;19.937;15.903;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0'
line =  '1.408;31.978;-3.772;0.217;30.077;-0.712;3.487;28.223;-0.125;5.525;31.454;-0.244;2.976;33.180;2.017;3.870;30.639;4.741'

DrawPDB(line)
ConstructPDB('DeNovo1.pdb')
#filename = 'DeNovo2.pdb'
#RotatePDB(filename)

#os.system('pymol DeNovo1.pdb')
#os.system('pymol DeNovo2.pdb')
#os.system('pymol DeNovo3.pdb')
#os.system('rm DeNovo1.pdb')
#os.system('rm DeNovo2.pdb')
#os.system('rm DeNovo3.pdb')
