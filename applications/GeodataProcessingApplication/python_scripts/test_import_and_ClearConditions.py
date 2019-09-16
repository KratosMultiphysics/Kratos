"""
# IT'S ONLY A TEST. SCRIPT NOT COMPLETE 
"""


import KratosMultiphysics as Kratos
import KratosMultiphysics.GeodataProcessingApplication as KratosGeo
import KratosMultiphysics.FluidDynamicsApplication

from geo_importer import GeoImporter
from geo_mesher import GeoMesher
from geo_preprocessor import GeoPreprocessor
from geo_building import GeoBuilding

import time
import os

start_time = time.time()
# num_test = "29_september_import_mdpa_and_ClearConditions"
num_test = "50_september_import_and_CleanConditions_cylinder"
print("\n\nTEST ", num_test, "\n\n")

# we create a new folders for this test
if not os.path.exists("cfd_data/test_{}".format(num_test)):
	os.mkdir("cfd_data/test_{}".format(num_test))
if not os.path.exists("cfd_data/test_{}/gid_file".format(num_test)):
	os.mkdir("cfd_data/test_{}/gid_file".format(num_test))
if not os.path.exists("cfd_data/test_{}/gid_file/Binary".format(num_test)):
	os.mkdir("cfd_data/test_{}/gid_file/Binary".format(num_test))
if not os.path.exists("cfd_data/test_{}/gid_file/Ascii".format(num_test)):
	os.mkdir("cfd_data/test_{}/gid_file/Ascii".format(num_test))
if not os.path.exists("cfd_data/test_{}/stl_file".format(num_test)):
	os.mkdir("cfd_data/test_{}/stl_file".format(num_test))
if not os.path.exists("cfd_data/test_{}/analysis_file".format(num_test)):
	os.mkdir("cfd_data/test_{}/analysis_file".format(num_test))
if not os.path.exists("cfd_data/test_{}/mdpa_file".format(num_test)):
	os.mkdir("cfd_data/test_{}/mdpa_file".format(num_test))


# import terrain
preprocessor = GeoPreprocessor()
importer = GeoImporter()
mesher = GeoMesher()
building = GeoBuilding()


# # import domain from mdpa file
# importer._InitializeModelPart("test_model")
# terrain_model_part = importer.ModelPart

# model_part_in = "data/mdpa_file/check/08_Box_buildings_distance_field_3"
# KratosMultiphysics.ModelPartIO(model_part_in).ReadModelPart(terrain_model_part)

# building.SetGeoModelPart(terrain_model_part)
# building.CreateGidControlOutput("cfd_data/test_{}/gid_file/Ascii/08_Box_buildings_distance_field_3".format(num_test), "GiD_PostAscii")
# building.CreateGidControlOutput("cfd_data/test_{}/gid_file/Binary/08_Box_buildings_distance_field_3".format(num_test), "GiD_PostBinary")

# building.HasModelPart = True
# building.HasDistanceField = True
# building.SubtractBuildingMOD(2.0, 100.0, 0.1, "exponential", "ISOSURFACE", "true")			# interpolation = exponential; disc_type = ISOSURFACE; remove_internal_regions=true
# building.CreateGidControlOutput("cfd_data/test_{}/gid_file/Ascii/09_Box_buildings_subtracted_3_before_CleanConditions".format(num_test), "GiD_PostAscii")
# building.CreateGidControlOutput("cfd_data/test_{}/gid_file/Binary/09_Box_buildings_subtracted_3_before_CleanConditions".format(num_test), "GiD_PostBinary")
# # writing file mdpa
# mdpa_out_name = "cfd_data/test_{}/mdpa_file/09_Box_buildings_subtracted_3_before_CleanConditions".format(num_test)
# building.WriteMdpaOutput(mdpa_out_name)


# main_model_part = building.GetGeoModelPart()

# # we delete the condition if at least one node it is not in main model part
# KratosGeo.CleaningUtilities(main_model_part).CleanConditions()

# building.SetGeoModelPart(main_model_part)

# # we write gid file
# building.CreateGidControlOutput("cfd_data/test_{}/gid_file/Ascii/10_Box_buildings_subtracted_3_after_CleanConditions".format(num_test), "GiD_PostAscii")
# building.CreateGidControlOutput("cfd_data/test_{}/gid_file/Binary/10_Box_buildings_subtracted_3_after_CleanConditions".format(num_test), "GiD_PostBinary")
# # we write mdpa file
# mdpa_out_name = "cfd_data/test_{}/mdpa_file/10_Box_buildings_subtracted_3_after_CleanConditions".format(num_test)
# building.WriteMdpaOutput(mdpa_out_name)



# import domain from mdpa file
importer._InitializeModelPart("test_model")
terrain_model_part = importer.ModelPart

model_part_in = "data/mdpa_file/domain_13_sept_2019/08_test_34_cylinder"
KratosMultiphysics.ModelPartIO(model_part_in).ReadModelPart(terrain_model_part)

building.SetGeoModelPart(terrain_model_part)
# writing file mdpa
mdpa_out_name = "cfd_data/test_{}/mdpa_file/08_Box_buildings_distance_field_3_new".format(num_test)
building.WriteMdpaOutput(mdpa_out_name)

#################################################
building = GeoBuilding()
building.SetGeoModelPart(terrain_model_part)
building.HasModelPart = True
building.HasDistanceField = True
# subtract buildings
# building.SubtractBuildingMOD(0.2, 10.0, 0.1, "exponential", "ISOSURFACE", "true")			# interpolation = exponential; disc_type = ISOSURFACE; remove_internal_regions=true
building.SubtractBuildingMOD(2.0, 100.0, 0.1, "exponential", "ISOSURFACE", "true")			# interpolation = exponential; disc_type = ISOSURFACE; remove_internal_regions=true
building.CreateGidControlOutput("cfd_data/test_{}/gid_file/Ascii/09_Box_buildings_subtracted_3_before_CleanConditions".format(num_test), "GiD_PostAscii")
building.CreateGidControlOutput("cfd_data/test_{}/gid_file/Binary/09_Box_buildings_subtracted_3_before_CleanConditions".format(num_test), "GiD_PostBinary")
# # writing file mdpa
# mdpa_out_name = "cfd_data/test_{}/mdpa_file/09_Box_buildings_subtracted_3_before_CleanConditions".format(num_test)
# building.WriteMdpaOutput(mdpa_out_name)



# #
# origin_model_part = building.GetGeoModelPart()

# # we save only SKIN_ISOSURFACE
# current_model = KratosMultiphysics.Model()

# # skin
# importer._InitializeModelPart("SKIN")
# SkinModelPart = importer.ModelPart

# skin_smp = origin_model_part.GetSubModelPart("SKIN_ISOSURFACE")

# for node in skin_smp.Nodes:
# 	SkinModelPart.CreateNewNode(node.Id, node.X, node.Y, node.Z)
# for cond in skin_smp.Conditions:
# 	nodes = cond.GetNodes()
# 	SkinModelPart.CreateNewElement("Element2D3N",  cond.Id, [nodes[0].Id, nodes[1].Id, nodes[2].Id], skin_smp.GetProperties()[0])

# building.SetGeoModelPart(SkinModelPart)
# building.CreateGidControlOutput("cfd_data/test_{}/gid_file/Ascii/99_Skin_model_part".format(num_test), "GiD_PostAscii")
# building.CreateGidControlOutput("cfd_data/test_{}/gid_file/Binary/99_Skin_model_part".format(num_test), "GiD_PostBinary")

# # bottom
# importer._InitializeModelPart("BOTTOM")
# BottomModelPart = importer.ModelPart

# bottom_smp = origin_model_part.GetSubModelPart("BottomModelPart")

# nodes_list = []
# for node in bottom_smp.Nodes:
# 	BottomModelPart.CreateNewNode(node.Id, node.X, node.Y, node.Z)
# 	nodes_list.append(node.Id)
# for cond in bottom_smp.Conditions:
# 	nodes = cond.GetNodes()
# 	for node in nodes:
# 		if (node.Id not in nodes_list):
# 			break
# 	else:
# 		BottomModelPart.CreateNewElement("Element2D3N",  cond.Id, [nodes[0].Id, nodes[1].Id, nodes[2].Id], skin_smp.GetProperties()[0])

# building.SetGeoModelPart(BottomModelPart)
# building.CreateGidControlOutput("cfd_data/test_{}/gid_file/Ascii/98_Bottom_model_part".format(num_test), "GiD_PostAscii")
# building.CreateGidControlOutput("cfd_data/test_{}/gid_file/Binary/98_Bottom_model_part".format(num_test), "GiD_PostBinary")

# input("PAUSE")
# #






# we delete the condition if at least one node it is not in main model part
stop_01 = time.time()
KratosGeo.CleaningUtilities(terrain_model_part).CleanConditions()
stop_02 = time.time()
print("\nCleanConditions done in {}\n".format(stop_02 - stop_01))


# # we delete the condition in the angles
# stop_01 = time.time()
# KratosGeo.CleaningUtilities(terrain_model_part).CleanConditionsAngles()
# stop_02 = time.time()
# print("\nCleanConditionsAngles done in {}\n".format(stop_02 - stop_01))



######################################################################
# it is like CleanConditionsAngles but in python

building.SetGeoModelPart(terrain_model_part)

stop_01 = time.time()

origin_model_part = building.GetGeoModelPart()
skin_smp = origin_model_part.GetSubModelPart("SKIN_ISOSURFACE")
bottom_smp = origin_model_part.GetSubModelPart("BottomModelPart")

print("*** Number nodes in SKIN_ISOSURFACE: ", skin_smp.NumberOfNodes())
print("*** Number nodes in BottomModelPart: ", bottom_smp.NumberOfNodes())
print("*** Number conditions in BottomModelPart: ", bottom_smp.NumberOfConditions())

skin_nodes = []
for node in skin_smp.Nodes:
	# skin_nodes.append([node.X, node.Y, node.Z])
	skin_nodes.append(node)
print("\tNumber nodes in skin_nodes: ", len(skin_nodes))

for cond in bottom_smp.Conditions:
	count = 0
	nodes = cond.GetNodes()
	for node in nodes:
		if node in skin_nodes:
			count += 1
	if (count == 3):
		print("\t*** the Condition \"{}\" will be delete!".format(cond.Id))
		cond.Set(KratosMultiphysics.TO_ERASE, True)
	
origin_model_part.RemoveConditionsFromAllLevels(KratosMultiphysics.TO_ERASE)
building.SetGeoModelPart(origin_model_part)

stop_02 = time.time()
print("*** Time: ", (stop_02 - stop_01))
######################################################################



building.SetGeoModelPart(terrain_model_part)
building.CreateGidControlOutput("cfd_data/test_{}/gid_file/Ascii/10_Box_buildings_subtracted_3_after_CleanConditions".format(num_test), "GiD_PostAscii")
building.CreateGidControlOutput("cfd_data/test_{}/gid_file/Binary/10_Box_buildings_subtracted_3_after_CleanConditions".format(num_test), "GiD_PostBinary")
# writing file mdpa
mdpa_out_name = "cfd_data/test_{}/mdpa_file/10_Box_buildings_subtracted_3_after_CleanConditions".format(num_test)
building.WriteMdpaOutput(mdpa_out_name)


print("*** Time: ", time.time() - start_time)
print("\nTEST ", num_test, "END\n\n")

