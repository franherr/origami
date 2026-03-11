import cmath
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection


#
# PROGRAM origami.py
#
# Version: February 27, 2026
#
# Has two modes, vector and raster
#
# Vector mode computes backward iterates of real line segment [0,1] under origami map 
# of Equilateral unit triangle to itself, depending on a parameter "tau"
# and user specified number "iterates" of iterates, and then plots the result
#
# Raster mode computes *forward* iterates of pixels and checks if they land in [0,1]
# depending on parameters "tau" and "iterates" (and number of pixels per line) and then
# plots the result
#
	

# operations on triangles
def is_in_edge(p,seg):
	l1 = length(seg)
	l2 = length([seg[0],p])
	l3 = length([p,seg[1]])
	if(l2+l3 - l1 < 0.01):
		return(1)
	else:
		return(0)

def is_in_triangle(p,X):
	# triangle X has shape (3,2)
	o = []
	M = np.zeros((2,2))
	for i in range(3):
		# is p above the ith oriented edge?
		# Doh! This assumes triangle is positively oriented!!!
		M[0] = X[(i+1)%3] - X[i]
		M[1] = p - X[i]
		if np.linalg.det(M)<0:
			o.append(0)
		else:
			o.append(1)
	if o[0] == o[1] and o[1] == o[2]:
		return(1)
	else:
		return(0)
	
def rotate_point_in_edge_of_triangle(p,X,i):
	midpoint = (X[(i+1)%3] + X[i])/2
	q = midpoint*2 - p
	return(q)
	
	
def rotate_into_triangle(p,X):
	# triangle X has shape (3,2)
	M = np.zeros((2,2))
	qq = np.zeros(2)
	for i in range(3):
	# is p above the ith oriented edge?
		M[0] = X[(i+1)%3] - X[i]
		M[1] = p - X[i]
		if np.linalg.det(M)<0:
			#rotate in ith oriented edge
			midpoint = (X[(i+1)%3] + X[i])/2
			qq = midpoint*2 - p
			break
	return(qq)
		
def triangle_map(X1,X2,p):	
	# computes image of p under affine map of R^2 taking triangle X1 to X2
	M1 = np.zeros((2,2))
	M1[0] = X1[1]-X1[0]
	M1[1] = X1[2]-X1[0]
	M1 = np.transpose(M1)
	#columns of M1 are edges of X1
	M2 = np.zeros((2,2))
	M2[0] = X2[1]-X2[0]
	M2[1] = X2[2]-X2[0]
	M2 = np.transpose(M2)
	#columns of M2 are edges of X2
	L = np.zeros((2,2))
	L = M2 @ (np.linalg.inv(M1))
	q = np.zeros(2)
	q = (L @ (p - X1[0])) + X2[0]
	#if q is out of range, rotate it back!
	if is_in_triangle(q,EQUI) == 0:
		q = rotate_into_triangle(q,EQUI)
	return(q)
	
def unrotated_triangle_map(X1,X2,p):	
	# computes image of p under affine map of R^2 taking triangle X1 to X2
	M1 = np.zeros((2,2))
	M1[0] = X1[1]-X1[0]
	M1[1] = X1[2]-X1[0]
	M1 = np.transpose(M1)
	#columns of M1 are edges of X1
	M2 = np.zeros((2,2))
	M2[0] = X2[1]-X2[0]
	M2[1] = X2[2]-X2[0]
	M2 = np.transpose(M2)
	#columns of M2 are edges of X2
	L = np.zeros((2,2))
	L = M2 @ (np.linalg.inv(M1))
	q = np.zeros(2)
	q = (L @ (p - X1[0])) + X2[0]
	#don't check out of range!
	return(q)
	
def triangle_segment_map(X1,X2,seg):	
	new_seg = np.zeros((2,2))
	new_seg[0]=unrotated_triangle_map(X1,X2,seg[0]) #unrotated 
	new_seg[1]=unrotated_triangle_map(X1,X2,seg[1])
	return(new_seg)
	
# operations on segments

def length(seg):
	da = np.zeros(2)
	da = seg[1]-seg[0]
	return(np.sqrt(da[0]*da[0] + da[1]*da[1]))

def edge(X,i):
	# returns edge i of triangle X
	new_seg = np.zeros((2,2))
	new_seg[0]=X[i]
	new_seg[1]=X[(i+1)%3]
	return(new_seg)

def rotate_segment_in_edge_of_triangle(seg,X,i):
	new_seg = np.zeros((2,2))
	new_seg[0]=rotate_point_in_edge_of_triangle(seg[0],X,i)
	new_seg[1]=rotate_point_in_edge_of_triangle(seg[1],X,i)
	return(new_seg)
	
def ccw(A,B,C):
    return (C[1]-A[1]) * (B[0]-A[0]) >= (B[1]-A[1]) * (C[0]-A[0]) 

def do_segs_intersect(seg1,seg2):	
	#returns true if line segments seg1 and seg2 intersect
    return ccw(seg1[0],seg2[0],seg2[1]) != ccw(seg1[1],seg2[0],seg2[1]) and ccw(seg1[0],seg1[1],seg2[0]) != ccw(seg1[0],seg1[1],seg2[1])

def return_segs_intersect(seg1,seg2): 	# a1 a2, b1 b2
	#returns intersection of these two lines; doesn't check if they intersect
	da = seg1[1]-seg1[0]
	db = seg2[1]-seg2[0]
	dp = seg1[0]-seg2[0]
	dap = np.zeros(2)
	dap[0] = -da[1]
	dap[1] = da[0]
	denom = np.dot(dap,db) # dap[0]*db[0] + dap[1]*db[1] 	# dot( dap, db)
	num = np.dot(dap,dp) # dap[0]*dp[0] + dap[1]*dp[1]		# dot( dap, dp )
	return ((num / (0.00001+denom.astype(float)))*db + seg2[0])

def return_seg_tri_intersect(seg,X):
	intersection_points = []
	pruned_intersection_points = []
	for i in range(2):
		# add points in interior of triangle
		if is_in_triangle(seg[i],X) == 1:
			intersection_points.append(seg[i])
	for i in range(3):
		# add points on edges of triangle
		for j in range(2):
			if is_in_edge(seg[j],edge(X,i)):
				intersection_points.append(seg[j])
		# add points of transverse intersection with edges of triangles				
		if do_segs_intersect(seg,edge(X,i)):	
			intersection_points.append(return_segs_intersect(seg,edge(X,i)))
#	if(len(intersection_points)==1):
#		print("found",len(intersection_points),"point of intersection")
#		print("intersection point: ",intersection_points)
#		print("seg: ",seg)
#		print("triangle: ",X)
	#prune repeat points in case of degeneracy, eg segment intersects vertex
	for point in intersection_points:
		if all(length([point,other_point])>min_distance for other_point in pruned_intersection_points):
			 pruned_intersection_points.append(point)
	return(pruned_intersection_points)

#main program
#global constants
min_distance = 0.0001
alpha = np.sqrt(3)/2
EQUI = np.zeros((3,2))
EQUI[0][0] = 0
EQUI[0][1] = 0
EQUI[1][0] = 1
EQUI[1][1] = 0
EQUI[2][0] = 1/2
EQUI[2][1] = np.sqrt(3)/2

#input global parameters tau and number of iterates
tau = float(input("Enter tau: "))
#tau = 0.1	# in general, a parameter
iterates = int(input("Enter iterates: "))

output_style = int(input("Output style (0=vector, 1=raster): "))
if output_style == 1:
	pixel_number = int(input("Number of pixels (1000 will take a long time!!): "))

#derive constants sigma, kappa, theta and mu
sigma = (1-tau)/2
kappa = np.sqrt(np.square((1-sigma)*alpha) + np.square(sigma - (1-sigma)/2))
theta = np.arcsin(alpha*tau/kappa)
mu = -np.tan(theta)/2
print("tau",tau,"mu",mu)

# compute triangles A[0]..A[6] 
# better way: find coordinates of 9 vertices (as arrays) and lock corners to these
# Each triangle A[i] is a 3x2 array of coordinates of its vertices; 
# A[i][j] is the jth vertex of the ith triangle, and A[i][j][0], A[i][j][1] are the x and y coordinates of that vertex
A = np.zeros((7,3,2))
A[0][0][0] = 0
A[0][0][1] = 0
A[0][1][0] = sigma
A[0][1][1] = 0
A[0][2][0] = sigma/2
A[0][2][1] = sigma*alpha

A[1][0][0] = sigma/2
A[1][0][1] = sigma*alpha
A[1][1][0] = sigma
A[1][1][1] = 0
A[1][2][0] = (1-sigma)/2
A[1][2][1] = (1-sigma)*alpha

A[2][0][0] = (1-sigma)/2
A[2][0][1] = (1-sigma)*alpha
A[2][1][0] = sigma
A[2][1][1] = 0
A[2][2][0] = 1-(sigma/2)
A[2][2][1] = sigma*alpha

A[3][0][0] = (1-sigma)/2
A[3][0][1] = (1-sigma)*alpha
A[3][1][0] = 1-(sigma/2)
A[3][1][1] = sigma*alpha
A[3][2][0] = (1+sigma)/2
A[3][2][1] = (1-sigma)*alpha

A[4][0][0] = (1-sigma)/2
A[4][0][1] = (1-sigma)*alpha
A[4][1][0] = (1+sigma)/2
A[4][1][1] = (1-sigma)*alpha
A[4][2][0] = 1/2
A[4][2][1] = alpha

A[5][0][0] = sigma
A[5][0][1] = 0
A[5][1][0] = 1-sigma
A[5][1][1] = 0
A[5][2][0] = 1-(sigma/2)
A[5][2][1] = sigma*alpha

A[6][0][0] = 1-sigma
A[6][0][1] = 0
A[6][1][0] = 1
A[6][1][1] = 0
A[6][2][0] = 1-(sigma/2)
A[6][2][1] = sigma*alpha

# compute their image triangles B[0]..B[6]
# better way: find coordinates of 9 vertices (as arrays) and lock corners to these
# auxiliary points
pp = np.zeros(2)
pp[0] = 7/8 - (alpha*mu)/2 - (alpha*(alpha-mu))/2
pp[1] = (alpha-mu)/4 - alpha*(1/4 + alpha*mu)

qq = np.zeros(2)
qq[0] = 3/8 + (alpha*mu)/2 + (alpha*(alpha+mu))/2
qq[1] = alpha - (alpha+mu)/4 - alpha*(1/4 - alpha*mu)

B = np.zeros((7,3,2))
B[0][0][0] = 0
B[0][0][1] = 0
B[0][1][0] = 1/2
B[0][1][1] = mu
B[0][2][0] = 1/4 - (alpha*mu)
B[0][2][1] = (alpha+mu)/2

B[1][0][0] = 1/4 - (alpha*mu)
B[1][0][1] = (alpha+mu)/2
B[1][1][0] = 1/2
B[1][1][1] = mu
B[1][2][0] = 1/4 + (alpha*mu)
B[1][2][1] = (alpha-mu)/2

B[2][0][0] = 1/4 + (alpha*mu)
B[2][0][1] = (alpha-mu)/2
B[2][1][0] = 1/2
B[2][1][1] = mu
B[2][2][0] = 3/4 - (alpha*mu)
B[2][2][1] = (alpha-mu)/2

B[3][0][0] = 1/4 + (alpha*mu)
B[3][0][1] = (alpha-mu)/2
B[3][1][0] = 3/4 - (alpha*mu)
B[3][1][1] = (alpha-mu)/2
B[3][2] = qq

B[4][0][0] = 1/4 + (alpha*mu)
B[4][0][1] = (alpha-mu)/2
B[4][1] = qq
B[4][2][0] = 1/2 
B[4][2][1] = alpha

B[5][0][0] = 1/2
B[5][0][1] = mu
B[5][1] = pp
B[5][2][0] = 3/4 - (alpha*mu)
B[5][2][1] = (alpha-mu)/2

B[6][0] = pp
B[6][1][0] = 1
B[6][1][1] = 0
B[6][2][0] = 3/4 - (alpha*mu)
B[6][2][1] = (alpha-mu)/2

if output_style == 0:
	#vector output; inverse map on segments

	#initialize segment list to a single segment from (0,0) to (1,0)
	segment_list = []
	new_segment_list = []
	seg = np.zeros((2,2))
	point0 = np.zeros(2)
	point1 = np.zeros(2)
	point1[0]=1
	seg[0]=point0
	seg[1]=point1
	segment_list.append(seg)

	#main dynamical iteration loop
	for k in range(iterates):
		for segment in segment_list:
			for i in range(7):
				for j in range(4):
					if j < 3:
						segment_2 = rotate_segment_in_edge_of_triangle(segment,EQUI,j)
					else:
						segment_2 = segment
					segment_3 = return_seg_tri_intersect(segment_2,B[i])	# try to intersect segment_3 with B[i]
					if len(segment_3) >= 2:
						segment_4 = triangle_segment_map(B[i],A[i],segment_3)
						new_segment_list.append(segment_4)
				#		print("intersected segment with triangle ",i," and I got one!")
					#else:
						#print("intersected segment with triangle ",i," and it was empty")
		segment_list = []
		for seg in new_segment_list:
			#prune list; only keep segments which are long enough
			if length(seg)>min_distance:
				segment_list.append(seg)
			else:
				print("found a shorty!")
		new_segment_list = []
		print("depth: ",k, "segments: ", len(segment_list))
	
	#plot the output using mathplotlib
	lc = LineCollection(segment_list)
	fig, ax = plt.subplots()
	ax.add_collection(lc)

	plt.show()

else:
	# raster output; forward map on pixels
	
	p_init = np.zeros(2)
	p = np.zeros(2)
	q = np.zeros(2)
	xlist = []
	ylist = []
	for i in range(3):
		xlist.append(EQUI[i][0])
		ylist.append(EQUI[i][1])

	for i in range(pixel_number):
		if i%10 == 0:
			print("computing row",i,"of",pixel_number)
		for j in range(pixel_number-i):
			p_init[0] = i/pixel_number + j/(2*pixel_number)
			p_init[1] = j*alpha/pixel_number
			p = p_init
			for l in range(iterates):
				for k in range(7):
					if is_in_triangle(p,A[k])==1:
						q = triangle_map(A[k],B[k],p)
						break
				p = q
			if p[1]<0.01 and p[1]>-0.01:
				xlist.append(p_init[0])
				ylist.append(p_init[1])
	#	xlist.append(p[0])
	#	ylist.append(p[1])
	plt.scatter(xlist, ylist, s=0.1)
	plt.show()
