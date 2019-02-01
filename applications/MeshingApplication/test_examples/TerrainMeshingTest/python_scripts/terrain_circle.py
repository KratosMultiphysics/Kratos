from scipy import spatial
import math
import numpy as np
import os
import sys

# path file
path_file = os.path.abspath(__file__)		# absolute path of the file
name_dir = os.path.dirname(path_file)		# directory of the file
parent_dir = os.path.dirname(name_dir)		# parent directory

stl_file_name = parent_dir + os.sep + "stl_file" + os.sep + "terrain_2.stl"


'''
terrain mesh
'''
# read terrain model to extract verices and elements information
with open (stl_file_name) as read_file:
	# vertices
	vertex = {}					# here there are all vertices. key = (x, y, z); value = ID_vertex
	ID_vertex = 1				# initializing the index

	ids_all = []				# here there are only id nodes that have Z != 0.0

	# elements
	ID_elem = 1
	elem_info = []				# here there are the information about elements

	for row in read_file.readlines():										# reading lines into an array
		row = row.split()													# splitting the row in seperated parts
																			# row[0] = "vertex"; row[1] = x coordinate; row[2] = y coordinate; row[3] = z coordinate
		if (row[0] == "outer"):
			vertices_to_element = []

		elif (row[0] == "vertex"):
			X_coord, Y_coord, Z_coord = [float(coord) for coord in row[1:]]

			if (X_coord, Y_coord, Z_coord) not in vertex:					# CASE: vertex does not yet exist
				vertex[X_coord, Y_coord, Z_coord] = ID_vertex				# creating a new vertex

				if Z_coord != 0.0:
					ids_all.append(ID_vertex) 								# id vertex with Z != 0.0

				vertices_to_element.append(ID_vertex)						# here there are the vertices that I need to create the element
				ID_vertex += 1

			else:															# CASE: vertex does already exist
				vertices_to_element.append(vertex[X_coord, Y_coord, Z_coord])

		elif (row[0] == "endloop"):
			elem_info.append([ID_elem, [vertices_to_element[0], vertices_to_element[1], vertices_to_element[2]]])
			ID_elem += 1

X = []
Y = []
Z = []
coord_3D = {}		# dictionary where key: id; value: [X, Y, Z]
coord_2D = []		# array with XY coordinates of vertices [X, Y]
for k, v in vertex.items():
	if (k[2] == 0.0):	# if (Z == 0.0)
		# fill just coord_3D and not coord_2D
		# # coord_3D.append([v, k[0], k[1], k[2]])
		continue
	X.append(k[0])
	Y.append(k[1])
	coord_2D.append([k[0], k[1]])
	coord_3D[v] = [k[0], k[1], k[2]]		# {id: [X, Y, Z]}
	Z.append(k[2])

# calculate minumin and maximum value of X and Y coordinate
x_min = min(X)
x_max = max(X)

y_min = min(Y)
y_max = max(Y)

# calculate the centre of circle
x_center = (x_min + x_max)/2
y_center = (y_min + y_max)/2

# calculate the radius considering the smaller side of the rectangle
r_boundary = min((x_max-x_min), (y_max-y_min))/2	# radius of the domain
r_ground = r_boundary * 2/3							# 2/3 of the radius of the domain
r_buildings = r_boundary / 3						# 1/3 of the radius of the domain

# buffering zone
delta = r_boundary/1000			# evaluate this value!!!

# calculate id nodes that are inside the r_ground
points_no_smooth = np.array(coord_2D)		# points that do not require smooth
tree = spatial.KDTree(points_no_smooth)
idx_inside = tree.query_ball_point([x_center, y_center], r = r_ground)		# node ids that are inside 2/3 of the circle

# here there are ids that must to be smooth
idx_to_smooth = list(set(ids_all) - set(idx_inside))

z_min = min(Z)

del_id = []						# list where there are the ids to be deleted because are outside of r_boundary
for index in idx_to_smooth:
	x_curr = coord_3D[index][0]	# current X coordinate
	y_curr = coord_3D[index][1]	# current Y coordinate
	dist = math.sqrt(((x_center-x_curr)**2)+((y_center-y_curr)**2))		# calculate the distance between current node and center of circle

	if (dist > r_boundary-delta):
		del_id.append(index)
		continue 							# here beta = 0 because the node is outside the circle boundary
	else:
		Z_beta = (-(dist-r_ground) / (r_boundary-r_ground))+1		# calculate beta

	coord_3D[index][2] = (coord_3D[index][2] - z_min) * Z_beta + z_min

# delete the nodes from coord_3D that are outside the r_boundary
for i in del_id:
	del coord_3D[i]


# create an array with 360 degree (from 0 to 360) -> 2*pi / 360 = pi /180
# theta = np.arange(0, 2*np.pi, np.pi/180)
theta = np.arange(0, 2*np.pi, np.pi/90)

# define the nodes on circle boundary
x_circle = []
y_circle = []
for th in theta:
	x_circle.append(r_boundary * np.cos(th) + x_center)	# move the node along X to consider that the circle have not the centre in (0.0, 0.0)
	y_circle.append(r_boundary * np.sin(th) + y_center)	# move the node along Y to consider that the circle have not the centre in (0.0, 0.0)


#######################################################################################################################
###													T R I A N G L E													###
#######################################################################################################################
import triangle

# fill x and y lists with 2D coordinate
x = []
y = []
for id, coords in coord_3D.items():
	x.append(coords[0])
	y.append(coords[1])

# add circle's nodes into x and y lists
for i in range(len(x_circle)):
	x.append(x_circle[i])
	y.append(y_circle[i])

# fill pts list with the couple x, y
pts = np.vstack((x, y)).T
A = dict(vertices = pts)
B = triangle.triangulate(A)		# B {vertices: [...], triangles: [...], vertex_markers: [...]}

# create a string to fil OBJ file
string_to_write = ""
for k, v in coord_3D.items():
	string_to_write += "v {} {} {}\n".format(v[0], v[1], v[2])
	# string_to_write += "v {} {} {}\n".format(v[0], v[2], v[1])

for i in range(len(x_circle)):
	string_to_write += "v {} {} {}\n".format(x_circle[i], y_circle[i], z_min)
	# string_to_write += "v {} {} {}\n".format(x_circle[i], z_min, y_circle[i])

for f in B["triangles"]:
	string_to_write += "f {} {} {}\n".format(f[0]+1, f[1]+1, f[2]+1)


'''
obj output
'''
obj_file_out = parent_dir + os.sep + "obj_file" + os.sep + "terrain_out.obj"
with open(obj_file_out, mode = "w", encoding = "utf-8") as fileWrite:
	fileWrite.write(string_to_write)
