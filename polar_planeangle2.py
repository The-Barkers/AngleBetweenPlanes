import os, sys
import math
import statistics
import random


#This takes in the following 2 files of tab separated polar coordinates of the surfaces of two planes, 
#and returns the angle between them in milli radians.
#Example line of input file: 
#RP_7    184.20054833    80.56027828     1388.29685373
#with arguments 0 being an unused point name string, and arguments 1, 2, 3 being 
#azimuthal angle (degrees), polar angle (degrees), and distance (mm)
#This requires python 3

#This can either have the filenames hard coded below, or entered as command line arguments as follows. 
#$ python3 planeangle.py data2.txt data1.txt

filename1 = "polar_data1.txt" #Plane1
filename2 = "polar_data2.txt" #Plane2
Ur_const= .00762 #mm
Ur_ppm= 2.5
Uphi= math.radians(1./3600.)#1 arc second
Utheta= math.radians(1./3600.)#1 arc second

def load_data(filename, delimiter=','):
	#returns a list of (x,y,z) tuples.
	data = []
	with open(filename,'r') as datfile:
		for line in datfile.readlines():
			words = line.strip().split(delimiter)
			if len(words)<4:
				print("skipping line ",line, "len ", len(words) )
				continue
			data.append( (float(words[1]), float(words[2]),float(words[3]) ))
	return data

if len(sys.argv) >= 3:
	filename1 = sys.argv[1]
	filename2 = sys.argv[2]
print(f"opening {filename1} and {filename2}")

convert = True #When true, read in xyz, then convert that to polar to test polar procedure.
data1 = load_data(filename1)
data2 = load_data(filename2)
#data1 = load_data(filename1,'\t') #come in as theta(degrees), phi(degrees), R
#data2 = load_data(filename2,'\t')


#############################################################################
class Cartesian_Sums:
	def __init__(self, xyz_data):
	    self.SX = 0.
	    self.SY = 0.
	    self.SZ = 0.
	    self.SX2 = 0.
	    self.SY2 = 0.
	    self.SZ2 = 0.
	    self.SXY = 0.
	    self.SYZ = 0.
	    self.SXZ = 0.
	    self.fill(xyz_data)
	    self.make_m()
	    self.make_mdet()
	    self.make_abc()

	def fill(self,xyz_data):
		for tup in xyz_data:
			self.SX  += tup[0]
			self.SX2 += tup[0]**2
			self.SY  += tup[1]
			self.SY2 += tup[1]**2
			self.SZ  += tup[2]
			self.SZ2 += tup[2]**2
			self.SXZ += tup[0]*tup[2]
			self.SXY += tup[0]*tup[1]
			self.SYZ += tup[1]*tup[2]
	def make_m(self):
	    self.m12 = self.SYZ*self.SXZ - self.SZ2*self.SXY
	    self.m13 = self.SXY*self.SYZ - self.SY2*self.SXZ
	    self.m23 = self.SXY*self.SXZ - self.SX2*self.SYZ
	    self.m11 = self.SY2*self.SZ2 - self.SYZ**2 
	    self.m22 = self.SX2*self.SZ2 - self.SXZ**2 
	    self.m33 = self.SX2*self.SY2 - self.SXY**2
	def make_mdet(self):
	    self.m_determinate = self.SX2*self.SY2*self.SZ2 + 2*self.SXY*self.SYZ*self.SXZ - self.SX2*self.SYZ**2 - self.SY2*self.SXZ**2 - self.SZ2*self.SXY**2 
	def make_abc(self):
	    self.a = -self.m11*self.SX - self.m12*self.SY - self.m13*self.SZ
	    self.b = -self.m12*self.SX - self.m22*self.SY - self.m23*self.SZ
	    self.c = -self.m13*self.SX - self.m23*self.SY - self.m33*self.SZ
	    self.d = 1.
	def get_abc(self, prime = False):
	    a = self.a
	    b = self.b
	    c = self.c
	    d = self.d
	    if prime:
	        m_determinate = self.m_determinate 
	        a *= m_determinate 
	        b *= m_determinate 
	        c *= m_determinate 
	        d = m_determinate 
	    return ( a, b, c, d)

def make_abc(data, prime = False):
	S=Cartesian_Sums(data)
	return S.get_abc(prime)

def make_alpha_beta(abc): #TODO: include in class
	#return (alpha, beta, gamma, delta) tuple from an (a, b, c, ~) tuple such as the output of make_abc
	delta = 1.0/math.sqrt(abc[0]**2 + abc[1]**2 + abc[2]**2)
	a = abc[0]*delta
	b = abc[1]*delta
	c = abc[2]*delta
	return (a, b, c, delta)

def make_Phi(abc1, abc2, degrees = False):
	#Given the parameters for 2 planes, returns the angle between them.
	#if degrees, returns results in degrees, else radians.
	delta1 = 1.0/math.sqrt(abc1[0]**2 + abc1[1]**2 + abc1[2]**2)
	delta2 = 1.0/math.sqrt(abc2[0]**2 + abc2[1]**2 + abc2[2]**2)
	arg = delta1*delta2*(abc1[0]*abc2[0] + abc1[1]*abc2[1] + abc1[2]*abc2[2])
	Phi = math.acos(arg)
	if Phi > 0.5*math.pi:
		Phi = math.pi - Phi
	if degrees:
		Phi = math.degrees(Phi)
	return Phi

def make_dropout_stdev(data1, data2, degrees = False):
	Phi_results = []
	for i in range(len(data1)):
		dta1 = data1.copy()
		del dta1[i]
		Phi_results.append( make_Phi( make_abc(dta1), make_abc(data2), degrees ) )
	for i in range(len(data2)):
		dta2 = data2.copy()
		del dta2[i]
		Phi_results.append( make_Phi( make_abc(data1), make_abc(dta2), degrees ) )
	f = math.sqrt(len(Phi_results) -1)
	return statistics.stdev(Phi_results)*f




#poldata:
#a list of tuples: (r, phi, theta)

#poldataU: 
#make this a list of tuples: (r, phi, theta, Ur, Uphi, Utheta)
#where the xyz are the xyz computed from the r, phi, theta

#enum for poldata and poldataU:
eR = 0
ephi = 1
etheta = 2
eUr = 3
eUphi = 4
eUtheta = 5

def make_cartesian_from_polar(polar_datum):
    #takes in a polar tuple (r, phi, theta_azimuth) and returns (x,y,z)
    r = polar_datum[eR]
    phi = polar_datum[ephi]
    theta = polar_datum[etheta]
    x = r*math.sin(phi)*math.cos(theta)
    y = r*math.sin(phi)*math.sin(theta)
    z = r*math.cos(phi)
    return (x, y, z)
	

def make_poldataU(poldata, Ur_const, Ur_ppm, Uphi, Utheta):
	#takes in a poldata list, returns a poldataU list and "data" like list of xyz tuples.
	poldataU = []
	xyz_data = [] 
	i=0
	for polar_datum in poldata:
		r = polar_datum[eR]
		poldataU.append( (r, polar_datum[ephi], polar_datum[etheta], Ur_const + Ur_ppm*r*0.000001, Uphi, Utheta) )
		xyz = make_cartesian_from_polar(polar_datum)
		xyz_data.append( xyz )
		if i==0:
			print(f"make_poldataU does r {polar_datum[eR]} phi {polar_datum[ephi]} radians theta {polar_datum[etheta]} radians -> x {xyz[0]} y {xyz[1]} z {xyz[2]}")
		i+=1
	return poldataU, xyz_data


def make_xyzdata(poldataU):
	return [make_cartesian_from_polar(polar_datum) for polar_datum in poldataU]
				
def make_var_sum(SJ, SJ_bar, poldataU, xyz_data, Phi, csc_term):
	#takes in Cartesian things SJ (current plane), SJ_bar (the other one), 
	#as well as the data, both in the form of polar + uncertainty poldataU, and cartesian xyz_data. 
	#and the co-secant term for this plane, and the Phi
	#then returns the sum of Phi variances resulting from the uncertainty in the measurements 
	#returns the sum of variances 
	varriance_sum = 0.
	abcJ = SJ.get_abc(prime = True)
	abcJbar = SJ_bar.get_abc(prime = True)
	for pol, xyz in zip(poldataU, xyz_data):
		(X, Y, Z) = xyz
		R = pol[eR]
		Sphi = math.sin(pol[ephi])
		Cphi = math.cos(pol[ephi])
		Stheta = math.sin(pol[etheta])
		Ctheta = math.cos(pol[etheta])
		#d(abc)/dR
		VdR = (-Sphi*Ctheta,   -Sphi*Stheta,   -Cphi)
		m11dR = 2.*(Y*Sphi*Stheta*SJ.SZ2 + Z*Cphi*SJ.SY2 - 2.*R*Sphi*Cphi*Stheta*SJ.SYZ)
		m22dR = 2.*(X*Sphi*Ctheta*SJ.SZ2 + Z*Cphi*SJ.SX2 - 2.*R*Sphi*Cphi*Ctheta*SJ.SXZ)
		m33dR = 2.*(X*Sphi*Ctheta*SJ.SY2 + Y*Sphi*Stheta*SJ.SX2 - 2.*R*(Sphi**2)*Stheta*Ctheta*SJ.SXY)
		m12dR = 2.*(R*Sphi*Cphi*Stheta*SJ.SXZ + R*Sphi*Cphi*Ctheta*SJ.SYZ - Z*Cphi*SJ.SXY - R*(Sphi**2)*Ctheta*Stheta*SJ.SZ2)
		m13dR = 2.*(R*(Sphi**2)*Ctheta*Stheta*SJ.SYZ + R*Sphi*Cphi*Stheta*SJ.SXY - Y*Sphi*Stheta*SJ.SXZ - R*Sphi*Cphi*Ctheta*SJ.SY2)
		m23dR = 2.*(R*(Sphi**2)*Ctheta*Stheta*SJ.SXZ + R*Sphi*Cphi*Ctheta*SJ.SXY - X*Sphi*Ctheta*SJ.SYZ - R*Sphi*Cphi*Stheta*SJ.SX2)
		dAdR = m11dR*(-SJ.SX) + m12dR*(-SJ.SY) + m13dR*(-SJ.SZ) + SJ.m11*VdR[0] + SJ.m12*VdR[1] + SJ.m13*VdR[2]
		dBdR = m12dR*(-SJ.SX) + m22dR*(-SJ.SY) + m23dR*(-SJ.SZ) + SJ.m12*VdR[0] + SJ.m22*VdR[1] + SJ.m23*VdR[2]
		dCdR = m13dR*(-SJ.SX) + m23dR*(-SJ.SY) + m33dR*(-SJ.SZ) + SJ.m13*VdR[0] + SJ.m23*VdR[1] + SJ.m33*VdR[2]
		dPhidR = SJ.cot_term*(dAdR + dBdR + dCdR) + csc_term*(abcJbar[0]*dAdR + abcJbar[1]*dBdR + abcJbar[2]*dCdR)
		varriance_sum += (dPhidR * pol[eUr])**2
		#d(abc)/dPhi
		Vdphi = (-R*Cphi*Ctheta,   -R*Cphi*Stheta,   -R*Sphi )
		m11dphi = 2.*Y*R*Cphi*Stheta*SJ.SZ2 - 2.*Z*R*Sphi*SJ.SY2 - 2.*R*R*math.cos(2.*pol[ephi])*Stheta*SJ.SYZ #TODO: factor 2
		m22dphi = 2.*X*R*Cphi*Ctheta*SJ.SZ2 - 2.*Z*R*Sphi*SJ.SX2 - 2.*R*R*math.cos(2.*pol[ephi])*Ctheta*SJ.SXZ
		m33dphi = 2.*X*R*Cphi*Ctheta*SJ.SY2 - 2.*Y*R*Cphi*Stheta*SJ.SX2 - 2.*R*R*math.sin(2.*pol[ephi])*Ctheta*Stheta*SJ.SXY
		m12dphi = R*R*math.cos(2.*pol[ephi])*Stheta*SJ.SXZ + R*R*math.cos(2.*pol[ephi])*Ctheta*SJ.SYZ + 2.*Z*R*Sphi*SJ.SXY - R*R*math.sin(2.*pol[ephi])*SJ.SZ2
		m13dphi = R*R*math.sin(2.*pol[ephi])*Ctheta*Stheta*SJ.SYZ + R*R*math.cos(2.*pol[ephi])*Stheta*SJ.SXY - 2.*Y*R*Cphi*Stheta*SJ.SXZ - R*R*math.cos(2.*pol[ephi])*Ctheta*SJ.SY2

		m23dphi = R*R*math.sin(2.*pol[ephi])*Ctheta*Stheta*SJ.SXZ + R*R*math.cos(2.*pol[ephi])*Ctheta*SJ.SXY - 2.*X*R*Cphi*Ctheta*SJ.SYZ - R*R*math.cos(2.*pol[ephi])*Stheta*SJ.SX2

		dAdphi = m11dphi*(-SJ.SX) + m12dphi*(-SJ.SY) + m13dphi*(-SJ.SZ) + SJ.m11*Vdphi[0] + SJ.m12*Vdphi[1] + SJ.m13*Vdphi[2]
		dBdphi = m12dphi*(-SJ.SX) + m22dphi*(-SJ.SY) + m23dphi*(-SJ.SZ) + SJ.m12*Vdphi[0] + SJ.m22*Vdphi[1] + SJ.m23*Vdphi[2]
		dCdphi = m13dphi*(-SJ.SX) + m23dphi*(-SJ.SY) + m33dphi*(-SJ.SZ) + SJ.m13*Vdphi[0] + SJ.m23*Vdphi[1] + SJ.m33*Vdphi[2]
		dPhidphi = SJ.cot_term*(dAdphi + dBdphi + dCdphi) + csc_term*(abcJbar[0]*dAdphi + abcJbar[1]*dBdphi + abcJbar[2]*dCdphi)
		varriance_sum += (dPhidphi * pol[eUphi])**2
		#d(abc)/dTheta
		Vdtheta = (Y, -X, 0)
		m11dtheta = 2.*X*Y*SJ.SZ2 - 2.*X*Z*SJ.SYZ
		m22dtheta = -2.*X*Y*SJ.SZ2 + 2.*Y*Z*SJ.SXZ
		m33dtheta = 2.*X*Y*(SJ.SX2 - SJ.SY2) - 2.*(X*X - Y*Y)*SJ.SXY
		m12dtheta = X*Z*SJ.SXZ - Y*Z*SJ.SYZ - (X*X - Y*Y)*SJ.SZ2
		m13dtheta = (X*X - Y*Y)*SJ.SYZ + X*Z*SJ.SXY - 2.*X*Y*SJ.SXZ + Y*Z*SJ.SY2
		m23dtheta = (X*X - Y*Y)*SJ.SXZ + Y*Z*SJ.SXY + 2.*X*Y*SJ.SYZ - X*Z*SJ.SX2
		dAdtheta = m11dtheta*(-SJ.SX) + m12dtheta*(-SJ.SY) + m13dtheta*(-SJ.SZ) + SJ.m11*Vdtheta[0] + SJ.m12*Vdtheta[1] + SJ.m13*Vdtheta[2]
		dBdtheta = m12dtheta*(-SJ.SX) + m22dtheta*(-SJ.SY) + m23dtheta*(-SJ.SZ) + SJ.m12*Vdtheta[0] + SJ.m22*Vdtheta[1] + SJ.m23*Vdtheta[2]
		dCdtheta = m13dtheta*(-SJ.SX) + m23dtheta*(-SJ.SY) + m33dtheta*(-SJ.SZ) + SJ.m13*Vdtheta[0] + SJ.m23*Vdtheta[1] + SJ.m33*Vdtheta[2]
		dPhidtheta = SJ.cot_term*(dAdtheta + dBdtheta + dCdtheta) + csc_term*(abcJbar[0]*dAdtheta + abcJbar[1]*dBdtheta + abcJbar[2]*dCdtheta)
		varriance_sum += (dPhidtheta * pol[eUtheta])**2
	return varriance_sum 

def make_measurement_uncert_polar(poldataU1, xyz_data1, poldataU2, xyz_data2, variance = False):
	#Takes in the raw data, 
	S1 = Cartesian_Sums(xyz_data1)
	S2 = Cartesian_Sums(xyz_data2)

	abc1 = S1.get_abc(prime = True)
	abc2 = S2.get_abc(prime = True)
	Phi = make_Phi(abc1, abc2)
	S1.n_mag2 = abc1[0]**2 + abc1[1]**2 + abc1[2]**2 
	S2.n_mag2 = abc2[0]**2 + abc2[1]**2 + abc2[2]**2
	S1.cot_term = math.cos(Phi)/(math.sin(Phi)*S1.n_mag2)
	S2.cot_term = math.cos(Phi)/(math.sin(Phi)*S2.n_mag2)
	csc_term = -1./(math.sqrt(S1.n_mag2)*math.sqrt(S2.n_mag2)*math.sin(Phi))
	S1.V = (-S1.SX, -S1.SY, -S1.SZ)
	S2.V = (-S2.SX, -S2.SY, -S2.SZ)
	
	varriance_sum = make_var_sum(S1, S2, poldataU1, xyz_data1, Phi, csc_term)
	varriance_sum+= make_var_sum(S2, S1, poldataU2, xyz_data2, Phi, csc_term)
	if not variance:
		return math.sqrt(varriance_sum)
	else:
		 return varriance_sum 

#############################################################################
#For Monte Carlo study
def perterb(poldataU):
	poldataUc = []
	for polar_datum in poldataU:
		tup = ( polar_datum[0]+ random.normalvariate(0, polar_datum[eUr]),\
			polar_datum[1]+ random.normalvariate(0, polar_datum[eUphi]),\
			polar_datum[2]+ random.normalvariate(0, polar_datum[eUtheta]),\
			polar_datum[3], polar_datum[4], polar_datum[5] ) 
		poldataUc.append(tup)
	return poldataUc 

def make_measurement_uncert_polar_MC(poldataU1, poldataU2, Nthrows):
	#equivolent to make_measurement_uncert_polar but done through MC perterbations on the data.
	varriance = 0.
	abc1_mean = make_abc(make_xyzdata(poldataU1))
	abc2_mean = make_abc(make_xyzdata(poldataU2))
	phi_mean = make_Phi(abc1_mean, abc2_mean)
	for i in range(Nthrows):
		poldataU1c = perterb(poldataU1)
		poldataU2c = perterb(poldataU2)
		abc1= make_abc(make_xyzdata(poldataU1c))
		abc2= make_abc(make_xyzdata(poldataU2c))
		phi= make_Phi(abc1, abc2)
		varriance += (phi - phi_mean)**2
	varriance /= Nthrows
	return math.sqrt(varriance)
		

#############################################################################
def printXYZ(data):
	for d in data:
		print(f"{d[0]}\t\t{d[1]}\t\t{d[2]}")

#      ____             _     
#     / __ )___  ____ _(_)___ 
#    / __  / _ \/ __ `/ / __ \
#   / /_/ /  __/ /_/ / / / / /
#  /_____/\___/\__, /_/_/ /_/ 
#             /____/          

#poldata:
#a list of tuples: (r, phi, theta)

#poldataU: 
#make this a list of tuples: (r, phi, theta, Ur, Uphi, Utheta)
#where the xyz are the xyz computed from the r, phi, theta

#Load the data
def make_poldata(data):
	#reformat input data (theta_degrees, phi_degrees, R) into poldata: (R, phi_radians, theta_radians)
	dat = []	
	for point in data:
		dat.append( (point[0], math.radians(point[2]), math.radians(point[1])) )
		#dat.append( (point[2], math.radians(point[1]), math.radians(point[0])) )
	return dat

if convert:
    poldata1 = make_poldata(data1)
    poldata2 = make_poldata(data2)
else:
    poldata1 = data1
    poldata2 = data2

poldataU1, xyz_data1 = make_poldataU(poldata1, Ur_const, Ur_ppm, Uphi, Utheta)
poldataU2, xyz_data2 = make_poldataU(poldata2, Ur_const, Ur_ppm, Uphi, Utheta)
abc1 = make_abc(xyz_data1)
abc2 = make_abc(xyz_data2)

varriance_measurement_MC = make_measurement_uncert_polar_MC(poldataU1, poldataU2, 1000)

#print("xyz1 ")
#printXYZ(xyz_data1)
#print("xyz2 ")
#printXYZ(xyz_data2)

#print("abc1", abc1)
#print("abc2", abc2)

#Uncomment these to see the normalized plane parameters
#print("abcd1:", make_alpha_beta(abc1))
#print("abcd2:", make_alpha_beta(abc2))

#MAIN COMPUTATION 
stdev_dropout = make_dropout_stdev(xyz_data1, xyz_data2)
varriance_measurement = make_measurement_uncert_polar(poldataU1, xyz_data1, poldataU2, xyz_data2, True)
uncertainty = math.sqrt(stdev_dropout**2 + varriance_measurement)
reco_phi = make_Phi(abc1, abc2)
print( f"{1000.*reco_phi} +- {1000.*uncertainty} milli-radians (measurement uncertainty {1000.*math.sqrt(varriance_measurement)}, drop-out {1000.*stdev_dropout})" )
print(f"Monte Carlo estimate of measurement uncertainty: {varriance_measurement_MC *1000} miliradians")
