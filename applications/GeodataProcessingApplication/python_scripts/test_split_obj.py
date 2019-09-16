"""
	edit:		03 September 2019 -> fill a ModelPart with nodes and elements
"""

import KratosMultiphysics as Kratos
import KratosMultiphysics.GeodataProcessingApplication as KratosGeo
import KratosMultiphysics.FluidDynamicsApplication

from geo_importer import GeoImporter
from geo_preprocessor import GeoPreprocessor

import time
import os

num_test = 2


# we create a new folders for this test
if not os.path.exists("data/test_split_{}".format(num_test)):
	os.mkdir("data/test_split_{}".format(num_test))
if not os.path.exists("data/test_split_{}/stl_file".format(num_test)):
	os.mkdir("data/test_split_{}/stl_file".format(num_test))


importer = GeoImporter()
preproc = GeoPreprocessor()

# IMPORT
# obj_in = "data/buildings/buildings_barcellona_1.obj"
# obj_in = "data/buildings/buildings_100.obj"
# obj_in = "data/buildings/buildings_num_9.obj"
# obj_in = "data/buildings/buildings_num_6.obj"
obj_in = "data/buildings/buildings_num_20.obj"

# we fill the node_map and elem_map
node_map, elem_map = importer.ObjToPyMap(obj_in)	# we fill also the ModelPart
copy_elem_map = elem_map.copy()

geometries = {}		# key: number of the geometry; value: list with element ids which belong at this n-th geometry
geom_id = 1			# geometry id
dict_visited = []	# list with the elements in the dictionary already visited

key_list = list(elem_map.keys())	# it is necessary because we delete elements in the dictionary; so we can not iterate on dictionary directly
for faces in key_list:
	# print("*** faces ", faces)
	if (faces in dict_visited):
		continue
	dict_visited.append(faces)

	next_faces = [faces]	# the next faces that will be processed

	node_visited = []	# list with nodes already visited in this geometry
	elem_visited = []	# list with elements already visited in this geometry
	while (next_faces):	# this loop ends when next_faces is empty
		current_face = next_faces[0]			# get the first element in next_faces up to the last element. In next operations we subtract elements in next_faces
		for node in elem_map[current_face]:
			elem_visited.append(current_face)	# we update elem_visited
			dict_visited.append(current_face)
			if (node in node_visited):
				continue						# we go on if the node is already visited
			node_visited.append(node)			# we update node_visited
			for id, elem in elem_map.items():
				if (node in elem):
					next_faces.append(id)
		next_faces = list(set(next_faces) - set(elem_visited))	# we remove duplicates and the elements that are already visited

		del(elem_map[current_face])			# we delete current_face in elem_map

	geometries[geom_id] = list(set(elem_visited))			# key: id geometry; value: id elements of this geometry
	geom_id += 1


elem_map = copy_elem_map		# restore the elem_map
# for geom_id, geom_faces in geometries.items():
# 	list_to_write = []
# 	for id_f in geom_faces:
# 		nodes = elem_map[id_f]
# 		for node_id in nodes:
# 			coord = node_map[node_id]

# 			list_to_write.append([coord[0], coord[1], coord[2]])
	
# 	# we write a .STL file
# 	stl_file_out = "data/test_split_{}/stl_file/building_{}.stl".format(num_test, geom_id)	# one file per building is created
# 	preproc.WriteSTL(stl_file_out, list_to_write)



# create 03/09/2019
# [ITA] viene creato un file obj in cui ogni volta che c'è un nuovo edificio viene scritto "o Building"
building_model_part = importer.GetGeoModelPart()
for geom_id, geom_faces in geometries.items():
	# print("geom_id ", geom_id)
	# print("geom_faces ", geom_faces)
	# print("***********************************************************")

	name_sub_model_building = "Building_{}".format(geom_id)
	current_sub_model_building = building_model_part.CreateSubModelPart(name_sub_model_building)

	for id_f in geom_faces:
		nodes = elem_map[id_f]
		current_sub_model_building.AddNodes(nodes)
		element_i = building_model_part.GetElement(id_f)
		current_sub_model_building.AddElement(element_i, 0)

	# list_to_write = []
	# for id_f in geom_faces:
	# 	nodes = elem_map[id_f]
	# 	for node_id in nodes:
	# 		coord = node_map[node_id]

	# 		list_to_write.append([coord[0], coord[1], coord[2]])

preproc.SetGeoModelPart(building_model_part)

# writing file mdpa
# mdpa_out_name = "data/mdpa_file/buildings_num_9"
# mdpa_out_name = "data/mdpa_file/buildings_num_6"
mdpa_out_name = "data/mdpa_file/buildings_num_20"
preproc.WriteMdpaOutput(mdpa_out_name)

# print(building_model_part)
