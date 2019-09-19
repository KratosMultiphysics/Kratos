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
num_test = "54_september_import_mdpa"
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
# if not os.path.exists("cfd_data/test_{}/analysis_file".format(num_test)):
# 	os.mkdir("cfd_data/test_{}/analysis_file".format(num_test))
if not os.path.exists("cfd_data/test_{}/mdpa_file".format(num_test)):
	os.mkdir("cfd_data/test_{}/mdpa_file".format(num_test))


# import terrain
preprocessor = GeoPreprocessor()
importer = GeoImporter()
mesher = GeoMesher()
extract_center = True


print("\n*** IMPORT TERRAIN FROM MDPA FILE ***\n")

# import domain from mdpa file
importer._InitializeModelPart("test_model")
terrain_model_part = importer.ModelPart
# model_part_in = "data/mdpa_file/05_Box_buildings_subtracted_1"
# model_part_in = "data/mdpa_file/09_Box_buildings_subtracted_3_before_CleanConditions"
model_part_in = "data/mdpa_file/check/11_Box_buildings_after_CleanConditionsAngles"
KratosMultiphysics.ModelPartIO(model_part_in).ReadModelPart(terrain_model_part)


print(terrain_model_part)


# sub_model_name = "SKIN_ISOSURFACE"

# for node in terrain_model_part.Nodes:
# 	node.Set(KratosMultiphysics.TO_ERASE, True)
# for node in terrain_model_part.GetSubModelPart(sub_model_name).Nodes:
# 	node.Set(KratosMultiphysics.TO_ERASE, False)
# terrain_model_part.RemoveNodesFromAllLevels(KratosMultiphysics.TO_ERASE)

# for cond in terrain_model_part.Conditions:
# 	cond.Set(KratosMultiphysics.TO_ERASE, True)
# for cond in terrain_model_part.GetSubModelPart(sub_model_name).Conditions:
# 	cond.Set(KratosMultiphysics.TO_ERASE, False)
# terrain_model_part.RemoveConditionsFromAllLevels(KratosMultiphysics.TO_ERASE)

# for elem in terrain_model_part.Elements:
# 	elem.Set(KratosMultiphysics.TO_ERASE, True)
# for elem in terrain_model_part.GetSubModelPart(sub_model_name).Elements:
# 	elem.Set(KratosMultiphysics.TO_ERASE, False)
# terrain_model_part.RemoveElementsFromAllLevels(KratosMultiphysics.TO_ERASE)





# terrain_model_part.RemoveSubModelPart("BottomModelPart")
# terrain_model_part.RemoveSubModelPart("Inlet")
# terrain_model_part.RemoveSubModelPart("Outlet")
# terrain_model_part.RemoveSubModelPart("Parts_Fluid")
# terrain_model_part.RemoveSubModelPart("TopModelPart")
# terrain_model_part.RemoveSubModelPart("SKIN_ISOSURFACE")


mesher.SetGeoModelPart(terrain_model_part)
mesher.CreateGidControlOutput("cfd_data/test_{}/gid_file/Ascii/11_Box_buildings_after_CleanConditionsAngles".format(num_test), "GiD_PostAscii")
mesher.CreateGidControlOutput("cfd_data/test_{}/gid_file/Binary/11_Box_buildings_after_CleanConditionsAngles".format(num_test), "GiD_PostBinary")
# writing file mdpa
mdpa_out_name = "cfd_data/test_{}/mdpa_file/11_Box_buildings_after_CleanConditionsAngles".format(num_test)
mesher.WriteMdpaOutput(mdpa_out_name)


input("PAUSE (press to continue)")

# if there is CenterCondition we create center_model_part
if (terrain_model_part.HasSubModelPart("CenterCondition")):
	print("\n*** CenterCondition ***\n")
	importer._InitializeModelPart("center_model_part")
	importer.HasModelPart = True
	center_model_part = importer.GetGeoModelPart()
	prop = center_model_part.Properties[0]
	sub_model_center = mesher.GetGeoModelPart().GetSubModelPart("CenterCondition")
	for node in sub_model_center.Nodes:
		center_model_part.CreateNewNode(node.Id, node.X, node.Y, node.Z)
	for cond in sub_model_center.Conditions:
		nodes = cond.GetNodes()
		# center_model_part.CreateNewCondition("SurfaceCondition3D3N", cond.Id, [nodes[0].Id, nodes[1].Id, nodes[2].Id], prop)
		# center_model_part.CreateNewElement("Element3D3N", cond.Id, [nodes[0].Id, nodes[1].Id, nodes[2].Id], prop)
		center_model_part.CreateNewElement("Element2D3N", cond.Id, [nodes[0].Id, nodes[1].Id, nodes[2].Id], prop)	# we need this to shift the buildings on terrain

# we clear the model part to import the domain after the MMG process
importer._InitializeModelPart("test_model_after_MMG")
terrain_model_part = importer.ModelPart
model_part_in = "data/mdpa_file/03_Mesh_cylinder_refinement_1"
KratosMultiphysics.ModelPartIO(model_part_in).ReadModelPart(terrain_model_part)

mesher.SetGeoModelPart(terrain_model_part)
mesher.CreateGidControlOutput("cfd_data/test_{}/gid_file/Ascii/03_Mesh_cylinder_refinement_1".format(num_test), "GiD_PostAscii")
mesher.CreateGidControlOutput("cfd_data/test_{}/gid_file/Binary/03_Mesh_cylinder_refinement_1".format(num_test), "GiD_PostBinary")





# for name_mdpa in range(1, 11):
# 	if (name_mdpa == 6) or (name_mdpa == 7):
# 		continue
# 	print("\n*** IMPORT \"{}\" FROM MDPA FILE ***\n".format(name_mdpa))

# 	# import domain from mdpa file
# 	importer._InitializeModelPart(str(name_mdpa))
# 	domain_model_part = importer.ModelPart
	
# 	model_part_in = "data/mdpa_file/check/{}".format(name_mdpa)
# 	KratosMultiphysics.ModelPartIO(model_part_in).ReadModelPart(domain_model_part)

# 	mesher.SetGeoModelPart(domain_model_part)
# 	mesher.CreateGidControlOutput("cfd_data/test_{}/gid_file/Ascii/{}".format(num_test, name_mdpa), "GiD_PostAscii")
# 	mesher.CreateGidControlOutput("cfd_data/test_{}/gid_file/Binary/{}".format(num_test, name_mdpa), "GiD_PostBinary")

# 	input("PAUSE (press to continue)")



# print("*** Time: ", time.time() - start_time)
# print("\nTEST ", num_test, "END\n\n")

