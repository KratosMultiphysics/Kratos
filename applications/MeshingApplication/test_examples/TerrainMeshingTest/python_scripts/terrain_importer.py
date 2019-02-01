
'''
TERRAIN MESH IMPORTER (file parser functions)
'''

########################################################################
# function imports the nodes from a *.stl file
# return: data container
#########################################################################

def stl_import( stl_file_name_input ):

	# read terrain model to extract vertices
	with open (stl_file_name_input) as read_file:

		vertex_dict = {}					# here there are all vertices. key = (x, y, z); value = ID_vertex
		vertex_id_list = []					# here there are only id nodes that have Z != 0.0
		id_vertex = 1						# initialization of the index counter

		coord_3D = {}		# dictionary where key: id; value: [X, Y, Z]
		coord_2D = []		# array with XY coordinates of vertices [X, Y]

		for row in read_file.readlines():

			row = row.split()

			if ( row[0] == "vertex" and all( isfloat(n) for n in row[1:] )):

				# row[0] = "vertex"; row[1] = x coordinate; row[2] = y coordinate; row[3] = z coordinate
				X_coord, Y_coord, Z_coord = [float(coord) for coord in row[1:]]

				if (X_coord, Y_coord, Z_coord) not in vertex_dict:
					vertex_dict[X_coord, Y_coord, Z_coord] = id_vertex

					if Z_coord != 0.0:
						vertex_id_list.append(id_vertex) # id vertex with Z != 0.0 (I exclude the nodes of the base)

					# here there are the vertices that I need to create the element
					id_vertex += 1

		for k, v in vertex_dict.items():
			if (k[2] == 0.0):	# if (Z == 0.0)
				continue
			coord_2D.append([k[0], k[1]])
			coord_3D[v] = [k[0], k[1], k[2]]		# {id: [X, Y, Z]}

	return ( vertex_id_list, vertex_dict, coord_2D, coord_3D )


########################################################################
# function imports the nodes from a *.xyz file
# return: data container
#########################################################################

def xyz_import( xyz_file_name_input ):

	# read terrain model to extract vertices
	with open (xyz_file_name_input) as read_file:

		vertex_dict = {}								# here there are all vertices. key = (x, y, z); value = ID_vertex
		vertex_id_list = []								# here there are only id nodes that have Z != 0.0
		id_vertex = 1									# initialization of the index counter

		coord_3D = {}		# dictionary where key: id; value: [X, Y, Z]
		coord_2D = []		# array with XY coordinates of vertices [X, Y]

		for row in read_file.readlines():

			row = row.split()

			if ( all( isfloat(n) for n in row) ):

				# row[0] = x coordinate; row[1] = y coordinate; row[2] = z coordinate
				X_coord, Y_coord, Z_coord = [float(coord) for coord in row[0:3]]
				vertex_dict[X_coord, Y_coord, Z_coord] = id_vertex

				if abs(Z_coord) > 0.0:
					vertex_id_list.append(id_vertex) 	# id vertex with Z != 0.0 (I exclude the nodes of the base)

				id_vertex += 1

		for k, v in vertex_dict.items():
			if (k[2] == 0.0):	# if (Z == 0.0)
				continue
			coord_2D.append([k[0], k[1]])
			coord_3D[v] = [k[0], k[1], k[2]]		# {id: [X, Y, Z]}

	return ( vertex_id_list, vertex_dict, coord_2D, coord_3D )


########################################################################
# auxiliary functions
#########################################################################

def isfloat(value):
  try:
    float( value )
    return True
  except ValueError:
    return False