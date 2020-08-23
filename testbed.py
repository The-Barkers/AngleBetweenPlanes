import math
import random
import os,sys

target_angle = 0.0001 #radians
interplane_distance = 0.5 #distance units
sample_square_size = 1 #distance units, side length of a square of measurements
n_points_per_side = 10

surface_roughness = 0.0001 #distance units
surface_warp_amplitude = 0.0001 #distance units
crinkliness = 1 #dimensionless integer > 0, represenitng the nubmer of wavelengths of warping within sample_square_size 

#which is to be defined in terms of the origin of plane1
tracker_location_XYZ = (1.,0., interplane_distance) #distance units
tracker_orientation_phi = 0.1 #radians range 0..pi.
tracker_orientation_theta = 0.1 #radians range 0..2pi, which direction to bend in phi. 0 is positive x axis.

#Generate ideal XY ponts on the X,Y plane
#Origin is on plane 1.
def generate_point_grid(size, n, Z):
	#generates a grid of points, n on a side, in a square of size 'size' 
	#on a horizontal plane at z=Z. 	
	points = []
	ds = size/(n-1)
	for i in range(n):
		x = -0.5*size + i*ds
		for j in range(n):
			y = -0.5*size + j*ds
			points.append((x,y,Z))
	return points
ideal_points_xyz1 = generate_point_grid(sample_square_size, n_points_per_side, 0)
ideal_points_xyz2 = generate_point_grid(sample_square_size, n_points_per_side, interplane_distance)

#introduce surface warping as vertical deviations sway
def warp(xyz, amplitude, patch_size, crinkliness  ): 
	points = []
	for c in xyz:
		dz = amplitude*(math.cos( math.pi * crinkliness * c[0]/patch_size) + math.cos( math.pi * crinkliness * c[1]/patch_size))
		points.append( (c[0],c[1],c[2]+dz) )
	return points 
warped_points_xyz1 = warp(ideal_points_xyz1, surface_warp_amplitude, sample_square_size , crinkliness )
warped_points_xyz2 = warp(ideal_points_xyz2, surface_warp_amplitude, sample_square_size , crinkliness )


#introduce surface roughness, as vertical jitter
def roughen(xyz, sr):
	points = []
	for c in xyz:
		points.append( (c[0], c[1], c[1]+random.normalvariate(0, sr) ) )
	return points 
roughed_points_xyz1 = roughen(warped_points_xyz1 , surface_roughness )
roughed_points_xyz2 = roughen(warped_points_xyz2 , surface_roughness )

#rotate points
def rotateXZ(xyz_data, Phi):
	#rotate coords by angle Phi in radians, pitching up 
	points = []
	for c in xyz_data:
		tup = (c[0]*math.cos(Phi) - c[2]*math.sin(Phi), c[1], c[0]*math.sin(Phi) + c[2]*math.cos(Phi))
		points.append( tup )
	return points 
rotated_points_xyz1 = roughed_points_xyz1
rotated_points_xyz2 = rotateXZ(roughed_points_xyz2, target_angle )


def write_to_file(data, filename, delimeter):
	with open(filename,'w') as f:
		i = 0
		for datum in data:
			f.write(f"p{i}{delimeter}{datum[0]}{delimeter}{datum[1]}{delimeter}{datum[2]}\n")
			i+=1
write_to_file(rotated_points_xyz1, "test_cartesian1.txt",',')
write_to_file(rotated_points_xyz2, "test_cartesian2.txt",',')

#reexpress all the coordinates in terms of the tracker coordinates by parrallel transport.
def translate_list(data, vec):
	points = []
	for c in data:
		tup = (c[0] - vec[0],c[1] - vec[1],c[2] - vec[2] )
		points.append( tup )
	return points 
	
translated_points_xyz1 = translate_list(rotated_points_xyz1, tracker_location_XYZ)
translated_points_xyz2 = translate_list(rotated_points_xyz2, tracker_location_XYZ)
#translated_points_xyz1 = (rotated_points_xyz1[0] - tracker_location_XYZ[0],rotated_points_xyz1[1] - tracker_location_XYZ[1],rotated_points_xyz1[2] - tracker_location_XYZ[2] )
#translated_points_xyz2 = (rotated_points_xyz2[0] - tracker_location_XYZ[0],rotated_points_xyz2[1] - tracker_location_XYZ[1],rotated_points_xyz2[2] - tracker_location_XYZ[2] )

#TODO rotate the coodinates 
def do_rotate(xyz, phi, theta ):
	return xyz #placeholder
tracker_cartesian1 = do_rotate(translated_points_xyz1, tracker_orientation_phi, tracker_orientation_theta )
tracker_cartesian2 = do_rotate(translated_points_xyz2, tracker_orientation_phi, tracker_orientation_theta )

#TODO reexpress the coordinates in terms of polar coordinates
def to_polar(xyz):
	points = []
	for c in xyz:
		r = math.sqrt(c[0]**2 + c[1]**2 + c[2]**2)
		phi = math.acos(c[2]/r)
		theta = math.atan(c[0]/c[1])
		tup = (r, phi, theta)
		points.append( tup )
	return points 
polar_ideal1 = to_polar(tracker_cartesian1)
polar_ideal2 = to_polar(tracker_cartesian2)

#introduce measurement error
def add_error(rphitheta, uRconst, uRppm, uPhi, uTheta):
	points = []
	for p in rphitheta:
		tup = (p[0] + random.normalvariate(0, uRconst + p[0]*0.000001*uRppm), p[1]+random.normalvariate(0, uPhi), p[2]+random.normalvariate(0, uTheta))
		points.append( tup ) 
	return points 
Ur_const= .00762 #mm
Ur_ppm= 2.5
Uphi= math.radians(1./3600.)#1 arc second
Utheta= math.radians(1./3600.)#1 arc seconds
real_polar1 = add_error(polar_ideal1, Ur_const, Ur_ppm, Uphi, Utheta)
real_polar2 = add_error(polar_ideal2, Ur_const, Ur_ppm, Uphi, Utheta)

def degreeify(polar):
	points = []
	for c in polar:
		tup = (c[0], math.degrees(c[1]), math.degrees(c[2]))
		points.append( tup )
	return points 

write_to_file(degreeify(real_polar1), "test_polar1.txt",',')
write_to_file(degreeify(real_polar2), "test_polar2.txt",',')

