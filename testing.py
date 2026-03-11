import cmath
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection

SIGFIGS = 4

tau = 1/3
alpha = np.sqrt(3)/2
EQUI = np.zeros((3,2))
EQUI[0][0] = 0
EQUI[0][1] = 0
EQUI[1][0] = 1
EQUI[1][1] = 0
EQUI[2][0] = 1/2
EQUI[2][1] = np.sqrt(3)/2

#derive constants sigma, kappa, theta and mu
sigma = (1-tau)/2
kappa = np.sqrt(np.square((1-sigma)*alpha) + np.square(sigma - (1-sigma)/2))
theta = np.arcsin(alpha*tau/kappa)
mu = -np.tan(theta)/2

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

def plotAtriangles(graph):
    for i in range(0,7):
        center = np.mean(A[i], axis=0)
        graph.text(center[0], center[1], str(i), fontsize=12, ha='center', va='center')
        for j in range(0,3):
            graph.plot([A[i][j][0], A[i][(j+1)%3][0]], [A[i][j][1], A[i][(j+1)%3][1]], 'k-')

def plotBtriangles(graph):
    for i in range(0,7):
        center = np.mean(B[i], axis=0)
        graph.text(center[0], center[1], str(i), fontsize=12, ha='center', va='center')
        for j in range(0,3):
            graph.plot([B[i][j][0], B[i][(j+1)%3][0]], [B[i][j][1], B[i][(j+1)%3][1]], 'r-')
            graph.plot(B[i][j][0], B[i][j][1], 'ro')

def plotSegment(graph, seg, vertices=False):
    graph.plot([seg[0][0], seg[1][0]], [seg[0][1], seg[1][1]], 'b-')
    if vertices:
        graph.plot(seg[0][0], seg[0][1], 'bo', markersize=3)
        graph.plot(seg[1][0], seg[1][1], 'bo', markersize=3)

def unrotated_triangle_map(X1,X2,p):	
	# computes image of p under affine map of R^2 taking triangle X1 to X2
    # confirmed that the orientation preserving/reversing is correct, and for the odd-numbered triangles,
    # the map reverses orientation.
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

#returns a boolean answer to whether a point q is in a triangle
#assumes that the triangle is not degenerate (three vertices are not colinear)
def pointInTriangle(q, tri):
    p0 = tri[0]
    p1 = tri[1]
    p2 = tri[2]

    if (p1[0] - p0[0] == 0):
        m01 = 'undefined'
    else:
        m01 = (p1[1] - p0[1])/(p1[0] - p0[0])
    if (p2[0] - p1[0] == 0):
        m12 = 'undefined'        
    else:       
        m12 = (p2[1] - p1[1])/(p2[0] - p1[0])
    if (p2[0] - p0[0] == 0):
        m20 = 'undefined'
    else:
        m20 = (p2[1] - p0[1])/(p2[0] - p0[0])

    if (m01 == 'undefined'):
        hs01 = q[0] - p0[0]
        side01 = p2[0] - p0[0]
    else:
        hs01 = m01*(q[0] - p0[0]) + p0[1] - q[1]
        side01 = m01*(p2[0] - p0[0]) + p0[1] - p2[1]

    if (m12 == 'undefined'):
        hs12 = q[0] - p1[0]
        side12 = p0[0] - p1[0]
    else:
        hs12 = m12*(q[0] - p1[0]) + p1[1] - q[1]
        side12 = m12*(p0[0] - p1[0]) + p1[1] - p0[1]
    
    if (m20 == 'undefined'):
        hs20 = q[0] - p2[0]
        side20 = p1[0] - p2[0]
    else: 
        hs20 = m20*(q[0] - p2[0]) + p2[1] - q[1]
        side20 = m20*(p1[0] - p2[0]) + p2[1] - p1[1]

    sameSide01 = (round(hs01*side01, SIGFIGS) >= 0)
    sameSide12 = (round(hs12*side12, SIGFIGS) >= 0)
    sameSide20 = (round(hs20*side20, SIGFIGS) >= 0)

    return (sameSide01 and sameSide12 and sameSide20)

#check if point (x,y) is in segment defined by endpoints seg[0] and seg[1]
#assumes that (x,y) is colinear with seg[0] and seg[1]
def pointInSegment(x,y,seg):
    if (x <= max(seg[0][0], seg[1][0]) and x >= min(seg[0][0], seg[1][0]) and y <= max(seg[0][1], seg[1][1]) and y >= min(seg[0][1], seg[1][1])):
        return True
    else:
        return False


#checks if seg and edge intersect
#if they are parallel or do not intersect, returns the point (10,10)
#otherwise, it returns their intersection point
def segmentsIntersect(seg, edge):
    #initialize variables
    segVert = False
    edgeVert = False
    slope_seg = 0
    slope_edge = 0
    intx = 10
    inty = 10

    #defines slope of segments and checks if they are vertical
    if (round(seg[0][0], SIGFIGS) == round(seg[1][0], SIGFIGS)):
        segVert = True
        slope_seg = 'undefined'
    else :
        slope_seg = (seg[1][1] - seg[0][1])/(seg[1][0] - seg[0][0])
    if (edge[0][0] == edge[1][0]):
        edgeVert = True
        slope_edge = 'undefined'
    else :
        slope_edge = (edge[1][1] - edge[0][1])/(edge[1][0] - edge[0][0])

    #returns if both segments vertical, computes intersection if only one of them is
    if (segVert and edgeVert):
        return [intx, inty]
    elif edgeVert :
        intx = edge[0][0]
        inty = slope_seg*(intx - seg[0][0]) + seg[0][1]
    elif segVert :
        intx = seg[0][0]
        inty = slope_edge*(intx - edge[0][0]) + edge[0][1]
    else:
        #computes intersection for nonparallel lines
        if (round(slope_seg, SIGFIGS) != round(slope_edge, SIGFIGS)):
            intx = (edge[0][1] - slope_edge*edge[0][0] - seg[0][1] + slope_seg*seg[0][0])/(slope_seg - slope_edge)
            inty = slope_seg*(intx - seg[0][0]) + seg[0][1]
        
    #checks if the intersection points are in both segments
    #if not, it sets them back to 10 and 10
    if (intx != 10 and inty != 10):
        if not(pointInSegment(intx, inty, seg) and pointInSegment(intx, inty, edge)):
            intx = 10
            inty = 10
    
    return [intx, inty]
            

#takes a segment and looks at the intersection with each triangle
#returns a list of 3-element lists consisting of a number for the triangle of intersection
#and the two endpoints of the corresponding pre-image segment
# e.g. [0, (0.5, 0), (0.25, 0.144)] represents an intersection with triangle B[0] with endpoints
# (0.5, 0) and (0.25, 0.144).
def divideSegment(seg):
    newSegments = []
  
    #adds new segments for each triangle
    for i in range(0,7):
       # print('checking triangle ' + str(i))
        intersections = set()

        #checks if the endpoints of the segment are inside the triangle
        if (pointInTriangle(seg[0], B[i])):
            intersections.add((seg[0][0], seg[0][1]))
        if (pointInTriangle(seg[1], B[i])):
            intersections.add((seg[1][0], seg[1][1]))

        #if both endpoints of the segment are in the triangle, it breaks the loop
        # and does not check intersections with edges
        if (len(intersections) > 1):
            newSegments.append([i] + list(intersections))
        else:
            for j in range(0,3):
                edge = [B[i][j], B[i][(j+1)%3]]
                intersection = segmentsIntersect(seg, edge)
                intx = intersection[0]
                inty = intersection[1]
                if (intx != 10 and inty != 10):
                    intersections.add((intx, inty))
            if (len(intersections) > 1):
                newSegments.append([i] + list(intersections))
            #  print('found segment intersection with triangle ' + str(i) + ': ' + str(intersections))
    return newSegments

def iteratePreimages(n):
    initialSegments = [[[0,0],[1,0]]]
    for i in range(0,n):
        preimageSegments =[]
        subdividedSegments = []
        for seg in initialSegments:
            divisions = divideSegment(seg)
            if (len(divisions) > 0):
                subdividedSegments += divisions
        for j in range(0, len(subdividedSegments)):
            p1 = unrotated_triangle_map(B[subdividedSegments[j][0]], A[subdividedSegments[j][0]], subdividedSegments[j][1])
            p2 = unrotated_triangle_map(B[subdividedSegments[j][0]], A[subdividedSegments[j][0]], subdividedSegments[j][2])
            preimageSegments.append([p1, p2])
        initialSegments = preimageSegments
        i += 1
    return initialSegments

def plotPreimages(n, graph):
    preimageSegments = iteratePreimages(n)
    for seg in preimageSegments:
        plotSegment(graph, seg)


fig, axes = plt.subplots(1, 2)

axes[0].set_aspect('equal')
axes[1].set_aspect('equal')
plotAtriangles(axes[0])
plotBtriangles(axes[1])
# plotAtriangles(axes[1])
# plotBtriangles(axes[2])

# print('triangle B[3] vertices: ' + str(B[3]))
# print('triangle B[4] vertices: ' + str(B[4]))

plotSegment(axes[1], [[0,0],[1,0]])
plotPreimages(10, axes[0])

# seg = [[0.25, 0.14433757], [0.33333333, 0.19245009]]
# plotSegment(axes[1], seg)
# divisions = divideSegment(seg)
# print(divisions)

# preimages = iteratePreimages(1)
# print(preimages[1][0])
# for i in range(0,6):
#     inside = pointInTriangle(preimages[1][0], B[i])
#     if inside :
#         print('point is inside triangle ' + str(i))
# # division = divideSegment(preimages[1])
# axes[1].plot([preimages[1][0][0], preimages[1][1][0]], [preimages[1][0][1], preimages[1][1][1]], 'g-')
# # print(division)

# seg = np.zeros((2,2))
# seg[1][0] = 1
# seg[1][1] = 1
# p1 = seg[0]
# p2 = seg[1]
# axes[1].plot([p1[0], p2[0]], [p1[1], p2[1]], 'b-')
# axes[1].text(p1[0], p1[1], str('p1'), fontsize=10)
# axes[1].text(p2[0], p2[1], str('p2'), fontsize=10)



#find all subdivided segments
# initialSegments = [[[0,0],[1,0]]]
# subdividedSegments = []
# for seg in initialSegments:
#     divisions = divideSegment(seg)
#     subdividedSegments += divisions
# print('number of subdivided segements: ' + str(len(subdividedSegments)))
# print(subdividedSegments)



# q1 = unrotated_triangle_map(A[1], B[1], p1)
# q2 = unrotated_triangle_map(A[1], B[1], p2)
# axes[1].plot([q1[0], q2[0]], [q1[1], q2[1]], 'b-')
# axes[1].text(q1[0], q1[1], str('q1'), fontsize=10)
# axes[1].text(q2[0], q2[1], str('q2'), fontsize=10)



plt.show()