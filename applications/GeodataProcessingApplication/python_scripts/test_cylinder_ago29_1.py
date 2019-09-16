""" cylindrical domain
	# create:	29 August 2019 (copy from test_cylinder_14_import_mdpa.py)

	# we can import the domain mdpa file (with refinement already done)
	# otherwise  we can import the STL file and run the MMG process
	# 
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
num_test = "august_03"
print("\n\nTEST ", num_test, "\n\n")

# we create a new folders for this test
if not os.path.exists("cfd_data/test_{}".format(num_test)):
	os.mkdir("cfd_data/test_{}".format(num_test))
if not os.path.exists("cfd_data/test_{}/gid_file".format(num_test)):
	os.mkdir("cfd_data/test_{}/gid_file".format(num_test))
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

import_terrain_mdpa = False
if import_terrain_mdpa:
	print("\n*** IMPORT TERRAIN FROM MDPA FILE ***\n")

	# import domain from mdpa file
	importer._InitializeModelPart("test_model")
	terrain_model_part = importer.ModelPart
	model_part_in = "data/mdpa_file/01_Mesh_cylinder"
	KratosMultiphysics.ModelPartIO(model_part_in).ReadModelPart(terrain_model_part)

	mesher.SetGeoModelPart(terrain_model_part)
	mesher.CreateGidControlOutput("cfd_data/test_{}/gid_file/01_Mesh_cylinder".format(num_test))

	# if there is CenterCondition we create center_model_part
	if (terrain_model_part.HasSubModelPart("CenterCondition")):
		print("*** CenterCondition ***")
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
	mesher.CreateGidControlOutput("cfd_data/test_{}/gid_file/03_Mesh_cylinder_refinement_1".format(num_test))

else:
	print("\n*** IMPORT TERRAIN FROM STL FILE ***\n")

	# import STL terrain and compute mesh circle
	stl_name = "terrain_barcelona_1"
	preprocessor.ReadSTL("data/terrain/{}.stl".format(stl_name))		# terrain_Barcelona_2

	X = []; Y = []
	for coord in preprocessor.point_list:
		X.append(coord[0])
		Y.append(coord[1])

	x_shift = max(X)/2
	y_shift = max(Y)/2
	preprocessor.Shift(-x_shift, -y_shift)

	# we extract the surface in the top of the terrain
	preprocessor.ExtractMountain(1.0)

	# we write the STL file shifted
	STL_name = "cfd_data/test_{}/stl_file/{}_shift.stl".format(num_test, stl_name)
	preprocessor.WriteSTL(STL_name)

	importer.StlImport(STL_name)
	terrain_model_part = importer.GetGeoModelPart()

	mesher.SetGeoModelPart(terrain_model_part)

	# we cut a circular portion of the terrain, we perform smoothing procedure and we compute the volume mesh
	height = 100.0		# height = 80.0
	# elem_list = mesher.MeshCircleWithTerrainPoints(height, 40, True)	# if the value is True we can save a list with the elements that are inside r_buildings
	mesher.MeshCircleWithTerrainPoints(height, 40, extract_center)	# if extract_center = True we create a sub model part with conditions that are inside r_buildings

	print(mesher.GetGeoModelPart())
	# input("CHECK \"CenterCondition\"")

	terrain_model_part = mesher.GetGeoModelPart()
	# writing GiD file
	mesher.CreateGidControlOutput("cfd_data/test_{}/gid_file/01_Mesh_cylinder".format(num_test))
	# writing file mdpa
	mdpa_out_name = "cfd_data/test_{}/mdpa_file/01_Mesh_cylinder".format(num_test)
	mesher.WriteMdpaOutput(mdpa_out_name)

	if extract_center:
		importer._InitializeModelPart("center_model_part")
		importer.HasModelPart = True	# IMPROVE IT! UGLY
		center_model_part = importer.GetGeoModelPart()
		prop = center_model_part.Properties[0]
		sub_model_center = mesher.GetGeoModelPart().GetSubModelPart("CenterCondition")
		for node in sub_model_center.Nodes:
			center_model_part.CreateNewNode(node.Id, node.X, node.Y, node.Z)
		for cond in sub_model_center.Conditions:
			nodes = cond.GetNodes()
			# center_model_part.CreateNewCondition("SurfaceCondition3D3N", cond.Id, [nodes[0].Id, nodes[1].Id, nodes[2].Id], prop)
			center_model_part.CreateNewElement("Element3D3N", cond.Id, [nodes[0].Id, nodes[1].Id, nodes[2].Id], prop)
	
	# we need to remove this sub model part to avoid problems in MMG process
	mesher.GetGeoModelPart().RemoveSubModelPart("CenterCondition")

	### MMG ###
	# 1st ground refinement
	print("\n\n***** 1st ground refinement *****\n")
	# distance field from ground
	mesher.ComputeDistanceFieldFromGround()
	mesher.CreateGidControlOutput("cfd_data/test_{}/gid_file/02_Mesh_cylinder_distance_field_1".format(num_test))
	# writing file mdpa
	mdpa_out_name = "cfd_data/test_{}/mdpa_file/02_Mesh_cylinder_distance_field_1".format(num_test)
	mesher.WriteMdpaOutput(mdpa_out_name)
	# refine mesh
	# mesher.RefineMesh_test(5.0, 10.0)
	mesher.RefineMesh_test(10.0, 100.0)
	mesher.CreateGidControlOutput("cfd_data/test_{}/gid_file/03_Mesh_cylinder_refinement_1".format(num_test))
	# writing file mdpa
	mdpa_out_name = "cfd_data/test_{}/mdpa_file/03_Mesh_cylinder_refinement_1".format(num_test)
	mesher.WriteMdpaOutput(mdpa_out_name)
	print("BOX REFINEMENT 1 DONE!")

# terrain volume mesh ready
main_model_part = mesher.GetGeoModelPart()


# only for check (delete it)
print(center_model_part)
mesher.SetGeoModelPart(center_model_part)
mesher.CreateGidControlOutput("cfd_data/test_{}/gid_file/center_model_part".format(num_test))
# input("\n03_Mesh_cylinder_refinement_1 DONE!\n")


# import buildings
# obj_file_in = "data/conference/Barcelona_2.obj"
# obj_file_in = "data/buildings/buildings.obj"
obj_file_in = "data/buildings/buildings_barcellona_1.obj"
importer.ObjImport(obj_file_in, "BuildingModelPart")
building_model_part = importer.GetGeoModelPart()

# print(building_model_part)

importer.CreateGidControlOutput("cfd_data/test_{}/gid_file/03_buildings".format(num_test))
# input("IMPORT BUILDINGS DONE!")

# STL_center_name = "cfd_data/test_{}/stl_file/terrain_center.stl".format(num_test)
# preprocessor.WriteSTL(STL_center_name, elem_list)
# importer.StlImport(STL_center_name)

# # terrain center
# terrain_model_part_center = importer.GetGeoModelPart()

building = GeoBuilding()
building.SetGeoModelPart(main_model_part)

if extract_center:
	building.ShiftBuildingOnTerrain(building_model_part, center_model_part)
	# building.DeleteBuildingsUnderValue(building_model_part, 10)	# ONLY FOR A CHECK
	building.CreateGidControlOutput("cfd_data/test_{}/gid_file/buildings_shifted".format(num_test))
	# writing file mdpa
	mdpa_out_name = "cfd_data/test_{}/mdpa_file/buildings_shifted".format(num_test)
	building.WriteMdpaOutput(mdpa_out_name)
input("PAUSE shift buildings")


##############################################################################
# building steps
print("\n\nSTART PROCESS BUILDING")
stop_begin_loop = time.time()

# building = GeoBuilding()
building.SetGeoModelPart(main_model_part)
building.ImportBuilding(building_model_part)

# 1st cut
print("\n\n***** STEP 1 *****\n")
# distance field from hull
building.ComputeDistanceFieldFromHull(False, 1e-7)
building.CreateGidControlOutput("cfd_data/test_{}/gid_file/04_Box_buildings_distance_field_1".format(num_test))
# writing file mdpa
mdpa_out_name = "cfd_data/test_{}/mdpa_file/04_Box_buildings_distance_field_1".format(num_test)
building.WriteMdpaOutput(mdpa_out_name)

input("\n04_Box_buildings_distance_field_1 DONE!")

# subtract buildings
building.SubtractBuildingMOD(0.5, 10.0, 0.1, "constant", "STANDARD", "false")			# interpolation = constant; disc_type = STANDARD; remove_internal_regions=false
building.CreateGidControlOutput("cfd_data/test_{}/gid_file/05_Box_buildings_subtracted_1".format(num_test))
# writing file mdpa
mdpa_out_name = "cfd_data/test_{}/mdpa_file/05_Box_buildings_subtracted_1".format(num_test)
building.WriteMdpaOutput(mdpa_out_name)

input("\n05_Box_buildings_subtracted_1 DONE!")

# # 2nd cut
# print("\n\n***** STEP 2 *****\n")
# # distance field from hull
# building.ComputeDistanceFieldFromHull(False, 1e-7)
# building.CreateGidControlOutput("cfd_data/test_{}/gid_file/06_Box_buildings_distance_field_2".format(num_test))
# # writing file mdpa
# mdpa_out_name = "cfd_data/test_{}/mdpa_file/06_Box_buildings_distance_field_2".format(num_test)
# building.WriteMdpaOutput(mdpa_out_name)
# # subtract buildings
# building.SubtractBuildingMOD(0.2, 10.0, 0.1, "constant", "STANDARD", "false")			# interpolation = linear; disc_type = STANDARD; remove_internal_regions=false
# building.CreateGidControlOutput("cfd_data/test_{}/gid_file/07_Box_buildings_subtracted_2".format(num_test))
# # writing file mdpa
# mdpa_out_name = "cfd_data/test_{}/mdpa_file/07_Box_buildings_subtracted_2".format(num_test)
# building.WriteMdpaOutput(mdpa_out_name)

# 3rd cut
print("\n\n***** STEP 3 *****\n")
# distance field from hull
building.ComputeDistanceFieldFromHull(False, 1e-7)
building.CreateGidControlOutput("cfd_data/test_{}/gid_file/08_Box_buildings_distance_field_3".format(num_test))
# writing file mdpa
mdpa_out_name = "cfd_data/test_{}/mdpa_file/08_Box_buildings_distance_field_3".format(num_test)
building.WriteMdpaOutput(mdpa_out_name)

input("\n08_Box_buildings_distance_field_3")

# subtract buildings
building.SubtractBuildingMOD(0.2, 10.0, 0.1, "exponential", "ISOSURFACE", "true")			# interpolation = exponential; disc_type = ISOSURFACE; remove_internal_regions=true
building.CreateGidControlOutput("cfd_data/test_{}/gid_file/09_Box_buildings_subtracted_3".format(num_test))
# writing file mdpa
mdpa_out_name = "cfd_data/test_{}/mdpa_file/09_Box_buildings_subtracted_3_before_CleanConditions".format(num_test)
building.WriteMdpaOutput(mdpa_out_name)

input("\n09_Box_buildings_subtracted_3")

# we update main_model_part
main_model_part = building.GetGeoModelPart()

print()
# we delete the condition if at least one node it is not in main model part
KratosGeo.CleaningUtilities(main_model_part).CleanConditions()

# we fill the bottom model part with the new conditions after MMG process
KratosGeo.CleaningUtilities(main_model_part).FillBottom()

# writing file mdpa
mdpa_out_name = "cfd_data/test_{}/mdpa_file/10_Box_buildings_subtracted_3_after_CleanConditions".format(num_test)
building.WriteMdpaOutput(mdpa_out_name)


print("*** Time: ", time.time() - start_time)
print("\nTEST ", num_test, "END\n\n")

