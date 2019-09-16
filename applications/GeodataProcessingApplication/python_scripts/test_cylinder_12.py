""" cylindrical domain
	# create:	12 July 2019
	# edit:		15 July 2019
	# """


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
num_test = "cylinder_08_with_elem_list"
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
preprocessor.ReadSTL("data/terrain/terrain_Barcelona_2.stl")

X = []; Y = []
for coord in preprocessor.point_list:
	X.append(coord[0])
	Y.append(coord[1])

x_shift = max(X)/2
y_shift = max(Y)/2
preprocessor.Shift(-x_shift, -y_shift)

# we extract the surface in the top of the terrain
preprocessor.ExtractMountain(1.0)

STL_name = "cfd_data/test_{}/stl_file/terrain_shift.stl".format(num_test)
preprocessor.WriteSTL(STL_name)

importer = GeoImporter()
importer.StlImport(STL_name)
terrain_model_part = importer.GetGeoModelPart()

mesher = GeoMesher()
mesher.SetGeoModelPart(terrain_model_part)

# we cut a circular portion of the terrain, we perform smoothing procedure and we compute the volume mesh
height = 100.0		# height = 80.0
elem_list = mesher.MeshCircleWithTerrainPoints(height, 40, True)	# if the value is True we can save a list with the elements that are inside r_buildings

##################################################################################################
STL_center_name = "cfd_data/test_{}/stl_file/terrain_center.stl".format(num_test)
preprocessor.WriteSTL(STL_center_name, elem_list)
importer.StlImport(STL_center_name)
importer.CreateGidControlOutput("cfd_data/test_{}/gid_file/terrain_center".format(num_test))
input(importer.ModelPart)
##################################################################################################

terrain_model_part = mesher.GetGeoModelPart()
mesher.CreateGidControlOutput("cfd_data/test_{}/gid_file/01_Mesh_cylinder".format(num_test))

input("01_Mesh_cylinder DONE!")

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

input("03_Mesh_cylinder_refinement_1 DONE!")

# terrain volume mesh ready
main_model_part = mesher.GetGeoModelPart()


# import buildings
# obj file
obj_file_in = "data/conference/Barcelona_2.obj"

importer.ObjImport(obj_file_in, "BuildingModelPart")
building_model_part = importer.GetGeoModelPart()
importer.CreateGidControlOutput("cfd_data/test_{}/gid_file/03_buildings".format(num_test))

STL_center_name = "cfd_data/test_{}/stl_file/terrain_center.stl".format(num_test)
preprocessor.WriteSTL(STL_center_name, elem_list)
importer.StlImport(STL_center_name)

# terrain center
terrain_model_part_center = importer.GetGeoModelPart()

building = GeoBuilding()
building.SetGeoModelPart(main_model_part)

building.ShiftBuildingOnTerrain(building_model_part, terrain_model_part_center)
building.DeleteBuildingsUnderValue(building_model_part, 10)
building.CreateGidControlOutput("cfd_data/test_{}/gid_file/buildings_shifted".format(num_test))
input("PAUSE shift buildings")


##############################################################################
# building steps
print("\n\nSTART PROCESS BUILDING")
stop_begin_loop = time.time()

# building = GeoBuilding()
# building.SetGeoModelPart(main_model_part)
building.ImportBuilding(building_model_part)

# 1st cut
print("\n\n***** STEP 1 *****\n")
# distance field from hull
building.ComputeDistanceFieldFromHull(False, 1e-7)
building.CreateGidControlOutput("cfd_data/test_{}/gid_file/04_Box_buildings_distance_field_1".format(num_test))
# writing file mdpa
mdpa_out_name = "cfd_data/test_{}/mdpa_file/04_Box_buildings_distance_field_1".format(num_test)
building.WriteMdpaOutput(mdpa_out_name)

input("04_Box_buildings_distance_field_1 DONE!")

# subtract buildings
building.SubtractBuildingMOD(0.5, 10.0, 0.1, "constant", "STANDARD", "false")			# interpolation = constant; disc_type = STANDARD; remove_internal_regions=false
building.CreateGidControlOutput("cfd_data/test_{}/gid_file/05_Box_buildings_subtracted_1".format(num_test))
# writing file mdpa
mdpa_out_name = "cfd_data/test_{}/mdpa_file/05_Box_buildings_subtracted_1".format(num_test)
building.WriteMdpaOutput(mdpa_out_name)

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
# subtract buildings
building.SubtractBuildingMOD(0.2, 10.0, 0.1, "exponential", "ISOSURFACE", "true")			# interpolation = exponential; disc_type = ISOSURFACE; remove_internal_regions=true
building.CreateGidControlOutput("cfd_data/test_{}/gid_file/09_Box_buildings_subtracted_3".format(num_test))
# writing file mdpa
mdpa_out_name = "cfd_data/test_{}/mdpa_file/09_Box_buildings_subtracted_3_before_CleanConditions".format(num_test)
building.WriteMdpaOutput(mdpa_out_name)

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
