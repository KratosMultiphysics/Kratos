""" cylindrical domain
	# create:	29 August 2019 (copy from test_cylinder_14_import_mdpa.py)
	# edit:		30 August 2019

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
num_test = "august_04_buildings"
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


# import buildings
importer = GeoImporter()
# obj_file_in = "data/conference/Barcelona_2.obj"
# obj_file_in = "data/buildings/buildings.obj"
obj_file_in = "data/buildings/buildings_barcellona_1.obj"
importer.ObjImport(obj_file_in, "BuildingModelPart")
building_model_part = importer.GetGeoModelPart()

print(building_model_part)

importer.CreateGidControlOutput("cfd_data/test_{}/gid_file/03_buildings".format(num_test))
input("IMPORT BUILDINGS DONE!")

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

