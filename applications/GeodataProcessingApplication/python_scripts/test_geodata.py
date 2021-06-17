import KratosMultiphysics as Kratos
import KratosMultiphysics.GeodataProcessingApplication as KratosGeo
import KratosMultiphysics.FluidDynamicsApplication

from geo_importer import GeoImporter
from geo_mesher import GeoMesher
from geo_preprocessor import GeoPreprocessor
from geo_building import GeoBuilding
from geo_model import GeoModel

import time
import os

start_time = time.time()
# the name of the test
num_test = "006"
print("\n\n[DEBUG PY] TEST ", num_test, "\n\n")

# we create a new folders for this test
if not os.path.exists("cfd_data"):
	os.mkdir("cfd_data")
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

#####################
### T E R R A I N ###
#####################
# import terrain
preprocessor = GeoPreprocessor()
importer = GeoImporter()
mesher = GeoMesher()
extract_center = True

# INPUT PARAMETER
# import_terrain_mdpa = False
import_terrain_mdpa = True

if import_terrain_mdpa:
	print("\n[DEBUG PY] IMPORT TERRAIN FROM MDPA FILE\n")

	# import domain from mdpa file
	importer._InitializeModelPart("test_model")
	terrain_model_part = importer.ModelPart
	model_part_folder = "data/mdpa_file/test_geodata/"
	# mdpa model with 1st volume mesh
	model_part_in = model_part_folder + "01_Mesh_cylinder"
	
	KratosMultiphysics.ModelPartIO(model_part_in).ReadModelPart(terrain_model_part)

	mesher.SetGeoModelPart(terrain_model_part)
	mesher.CreateGidControlOutput("cfd_data/test_{}/gid_file/Ascii/01_Mesh_cylinder".format(num_test), "GiD_PostAscii")
	mesher.CreateGidControlOutput("cfd_data/test_{}/gid_file/Binary/01_Mesh_cylinder".format(num_test), "GiD_PostBinary")

	# if there is CenterCondition we create center_model_part
	if (terrain_model_part.HasSubModelPart("CenterCondition")):
		KratosMultiphysics.Logger.PrintWarning("test_{}".format(num_test), "sub model part \"CenterCondition\" found.")
		importer._InitializeModelPart("center_model_part")
		importer.HasModelPart = True
		center_model_part = importer.GetGeoModelPart()
		prop = center_model_part.Properties[0]
		sub_model_center = mesher.GetGeoModelPart().GetSubModelPart("CenterCondition")
		for node in sub_model_center.Nodes:
			center_model_part.CreateNewNode(node.Id, node.X, node.Y, node.Z)
		for cond in sub_model_center.Conditions:
			nodes = cond.GetNodes()
			center_model_part.CreateNewElement("Element2D3N", cond.Id, [nodes[0].Id, nodes[1].Id, nodes[2].Id], prop)	# we need this to shift the buildings on terrain
	
	# we clear the model part to import the domain after the MMG process
	importer._InitializeModelPart("test_model_after_MMG")
	terrain_model_part = importer.ModelPart
	# mdpa model with 1st refinement
	model_part_in = model_part_folder + "03_Mesh_cylinder_refinement_1"

	KratosMultiphysics.ModelPartIO(model_part_in).ReadModelPart(terrain_model_part)

	mesher.SetGeoModelPart(terrain_model_part)
	mesher.CreateGidControlOutput("cfd_data/test_{}/gid_file/Ascii/03_Mesh_cylinder_refinement_1".format(num_test), "GiD_PostAscii")
	mesher.CreateGidControlOutput("cfd_data/test_{}/gid_file/Binary/03_Mesh_cylinder_refinement_1".format(num_test), "GiD_PostBinary")

else:
	print("\n[DEBUG PY] IMPORT TERRAIN FROM STL FILE\n")

	# import STL terrain and compute mesh circle
	stl_name = "terrain_test"
	preprocessor.ReadSTL("data/terrain/{}.stl".format(stl_name))

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
	height = 100.0
	# num_sector = 12		# now in json file
	num_sectors = preprocessor.file_param["problem_data"]["num_sectors"].GetInt()
	print("\n\n[DEBUG PY] number of sectors: {}\n".format(num_sectors))
	
	mesher.MeshCircleWithTerrainPoints_old(height, 40, num_sectors, extract_center)	# if extract_center=True we create a sub model part with conditions that are inside r_buildings

	# we populate the sub model part sectors
	mesher.MeshSectors(num_sectors)

	# writing GiD file
	mesher.CreateGidControlOutput("cfd_data/test_{}/gid_file/Ascii/01_Mesh_cylinder".format(num_test), "GiD_PostAscii")
	mesher.CreateGidControlOutput("cfd_data/test_{}/gid_file/Binary/01_Mesh_cylinder".format(num_test), "GiD_PostBinary")
	
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
			center_model_part.CreateNewElement("Element2D3N", cond.Id, [nodes[0].Id, nodes[1].Id, nodes[2].Id], prop)
		
	# we need to remove this sub model part to avoid problems in MMG process
	mesher.GetGeoModelPart().RemoveSubModelPart("CenterCondition")

	### MMG ###
	# 1st ground refinement
	print("\n\n[DEBUG PY] 1st ground refinement\n")
	# distance field from ground
	mesher.ComputeDistanceFieldFromGround()
	mesher.CreateGidControlOutput("cfd_data/test_{}/gid_file/Ascii/02_Mesh_cylinder_distance_field_1".format(num_test), "GiD_PostAscii")
	mesher.CreateGidControlOutput("cfd_data/test_{}/gid_file/Binary/02_Mesh_cylinder_distance_field_1".format(num_test), "GiD_PostBinary")
	# writing file mdpa
	mdpa_out_name = "cfd_data/test_{}/mdpa_file/02_Mesh_cylinder_distance_field_1".format(num_test)
	mesher.WriteMdpaOutput(mdpa_out_name)
	# refine mesh
	mesher.RefineMesh_test(10.0, 100.0)
	mesher.CreateGidControlOutput("cfd_data/test_{}/gid_file/Ascii/03_Mesh_cylinder_refinement_1".format(num_test), "GiD_PostAscii")
	mesher.CreateGidControlOutput("cfd_data/test_{}/gid_file/Binary/03_Mesh_cylinder_refinement_1".format(num_test), "GiD_PostBinary")
	# writing file mdpa
	mdpa_out_name = "cfd_data/test_{}/mdpa_file/03_Mesh_cylinder_refinement_1".format(num_test)
	mesher.WriteMdpaOutput(mdpa_out_name)
	print("[DEBUG PY] BOX REFINEMENT 1 DONE!")


print("\n[DEBUG PY] TERRAIN DONE!\n")

# terrain volume mesh ready
main_model_part = mesher.GetGeoModelPart()

print("[DEBUG PY]\n", center_model_part)
mesher.SetGeoModelPart(center_model_part)
mesher.CreateGidControlOutput("cfd_data/test_{}/gid_file/Ascii/center_model_part".format(num_test), "GiD_PostAscii")
mesher.CreateGidControlOutput("cfd_data/test_{}/gid_file/Binary/center_model_part".format(num_test), "GiD_PostBinary")


#######################
### B U I L D I N G ###
#######################
building = GeoBuilding()
building.SetGeoModelPart(main_model_part)

# import buildings
import_building_mdpa = False
if import_building_mdpa:
	print("\n[DEBUG PY] IMPORT BUILDINGS FROM MDPA FILE\n")
	# import domain from mdpa file
	model_part_in = "/data/mdpa_file/mdpa_buildings"
	building.ImportBuildingHullMDPA(model_part_in)
	building_model_part = building.GetBuildingModelPart()

else:
	print("\n[DEBUG PY] IMPORT BUILDINGS FROM OBJ FILE\n")
	# obj file name
	# obj_file_in = "data/buildings/buildings_barcellona_test.obj"
	obj_file_in = "data/buildings/buildings_num_20.obj"
	# importer.ObjImportBuildings(obj_file_in, "BuildingModelPart", True)		# (obj_file_name_input, name_model_part, change_coord)
	importer.ObjToSplit(obj_file_in)
	# importer.ObjImportBuildings(obj_file_in, "BuildingModelPart", False)		# (obj_file_name_input, name_model_part, change_coord)
	building_model_part = importer.GetGeoModelPart()
print("[DEBUG PY]\n", building_model_part)

# we write GiD file
mesher.SetGeoModelPart(building_model_part)
mesher.CreateGidControlOutput("cfd_data/test_{}/gid_file/Ascii/buildings".format(num_test), "GiD_PostAscii")
mesher.CreateGidControlOutput("cfd_data/test_{}/gid_file/Binary/buildings".format(num_test), "GiD_PostBinary")

# writing file mdpa
mdpa_out_name = "cfd_data/test_{}/mdpa_file/buildings".format(num_test)
mesher.WriteMdpaOutput(mdpa_out_name)

print("\n[DEBUG PY] IMPORT BUILDINGS DONE!")

if extract_center:
	building.SetGeoModelPart(building_model_part)
	building.ShiftBuildingOnTerrain(center_model_part)
	print("\n[DEBUG PY] ShiftBuildingOnTerrain DONE!")

	# delete building placed out of center_model_part. If a building is out of range of the center_model_part is placed on z = 0
	building.DeleteBuildingsUnderValue(z_value=10)
	building.CreateGidControlOutput("cfd_data/test_{}/gid_file/Ascii/buildings_shifted".format(num_test), "GiD_PostAscii")
	building.CreateGidControlOutput("cfd_data/test_{}/gid_file/Binary/buildings_shifted".format(num_test), "GiD_PostBinary")
	# writing file mdpa
	mdpa_out_name = "cfd_data/test_{}/mdpa_file/buildings_shifted".format(num_test)
	building.WriteMdpaOutput(mdpa_out_name)

# building steps
print("\n\n[DEBUG PY] START PROCESS BUILDING")
stop_begin_loop = time.time()

# building = GeoBuilding()
building.SetGeoModelPart(main_model_part)		# it is a duplicate. CHECK IT
building.ImportBuilding(building_model_part)

# # 1st cut
# print("\n\n[DEBUG PY] STEP 1\n")
# # distance field from hull
# building.ComputeDistanceFieldFromHull(False, 1e-7)
# building.CreateGidControlOutput("cfd_data/test_{}/gid_file/Ascii/04_Box_buildings_distance_field_1".format(num_test), "GiD_PostAscii")
# building.CreateGidControlOutput("cfd_data/test_{}/gid_file/Binary/04_Box_buildings_distance_field_1".format(num_test), "GiD_PostBinary")
# # writing file mdpa
# mdpa_out_name = "cfd_data/test_{}/mdpa_file/04_Box_buildings_distance_field_1".format(num_test)
# building.WriteMdpaOutput(mdpa_out_name)
# # subtract buildings
# building.SubtractBuildingMOD(5.0, 100.0, 0.1, "Linear", "STANDARD", "false")			# interpolation = constant; disc_type = STANDARD; remove_internal_regions=false
# building.CreateGidControlOutput("cfd_data/test_{}/gid_file/Ascii/05_Box_buildings_subtracted_1".format(num_test), "GiD_PostAscii")
# building.CreateGidControlOutput("cfd_data/test_{}/gid_file/Binary/05_Box_buildings_subtracted_1".format(num_test), "GiD_PostBinary")
# # writing file mdpa
# mdpa_out_name = "cfd_data/test_{}/mdpa_file/05_Box_buildings_subtracted_1".format(num_test)
# building.WriteMdpaOutput(mdpa_out_name)

# # 2nd cut
# print("\n\n[DEBUG PY] STEP 2\n")
# # distance field from hull
# building.ComputeDistanceFieldFromHull(False, 1e-7)
# building.CreateGidControlOutput("cfd_data/test_{}/gid_file/Ascii/06_Box_buildings_distance_field_2".format(num_test), "GiD_PostAscii")
# building.CreateGidControlOutput("cfd_data/test_{}/gid_file/Binary/06_Box_buildings_distance_field_2".format(num_test), "GiD_PostBinary")
# # writing file mdpa
# mdpa_out_name = "cfd_data/test_{}/mdpa_file/06_Box_buildings_distance_field_2".format(num_test)
# building.WriteMdpaOutput(mdpa_out_name)
# # subtract buildings
# building.SubtractBuildingMOD(3.0, 10.0, 0.1, "Constant", "STANDARD", "false")			# interpolation = linear; disc_type = STANDARD; remove_internal_regions=false
# building.CreateGidControlOutput("cfd_data/test_{}/gid_file/Ascii/07_Box_buildings_subtracted_2".format(num_test), "GiD_PostAscii")
# building.CreateGidControlOutput("cfd_data/test_{}/gid_file/Binary/07_Box_buildings_subtracted_2".format(num_test), "GiD_PostBinary")
# # writing file mdpa
# mdpa_out_name = "cfd_data/test_{}/mdpa_file/07_Box_buildings_subtracted_2".format(num_test)
# building.WriteMdpaOutput(mdpa_out_name)

# 3rd cut
print("\n\n[DEBUG PY] STEP 3\n")
# distance field from hull
building.ComputeDistanceFieldFromHull(False, 1e-7)
print("\n\n[DEBUG PY] ComputeDistanceFieldFromHull\n")
building.CreateGidControlOutput("cfd_data/test_{}/gid_file/Ascii/08_Box_buildings_distance_field_3".format(num_test), "GiD_PostAscii")
building.CreateGidControlOutput("cfd_data/test_{}/gid_file/Binary/08_Box_buildings_distance_field_3".format(num_test), "GiD_PostBinary")
# writing file mdpa
mdpa_out_name = "cfd_data/test_{}/mdpa_file/08_Box_buildings_distance_field_3".format(num_test)
building.WriteMdpaOutput(mdpa_out_name)
# subtract buildings
building.SubtractBuildingMOD(2.0, 100.0, 0.1, "exponential", "ISOSURFACE", "true")			# interpolation = exponential; disc_type = ISOSURFACE; remove_internal_regions=true
print("\n\n[DEBUG PY] SubtractBuildingMOD\n")
building.CreateGidControlOutput("cfd_data/test_{}/gid_file/Ascii/09_Box_buildings_subtracted_3_before_CleanConditions".format(num_test), "GiD_PostAscii")
building.CreateGidControlOutput("cfd_data/test_{}/gid_file/Binary/09_Box_buildings_subtracted_3_before_CleanConditions".format(num_test), "GiD_PostBinary")

# we update main_model_part
main_model_part = building.GetGeoModelPart()

stop_01 = time.time()
# we delete the condition if at least one node it is not in main model part
KratosGeo.CleaningUtilities(main_model_part).CleanConditions()
stop_02 = time.time()
print("[DEBUG PY] CleanConditions done in ", (stop_02 - stop_01))

# we set the DENSITY and DYNAMIC_VISCOSITY values
prop = main_model_part.GetProperties()[0]
prop.SetValue(KratosMultiphysics.DENSITY, 1)
prop.SetValue(KratosMultiphysics.DYNAMIC_VISCOSITY, 0.002)
building.SetGeoModelPart(main_model_part)

# we write gid file
building.CreateGidControlOutput("cfd_data/test_{}/gid_file/Ascii/10_Box_buildings_subtracted_3_after_CleanConditions".format(num_test), "GiD_PostAscii")
building.CreateGidControlOutput("cfd_data/test_{}/gid_file/Binary/10_Box_buildings_subtracted_3_after_CleanConditions".format(num_test), "GiD_PostBinary")
# we write mdpa file
mdpa_out_name = "cfd_data/test_{}/mdpa_file/10_Box_buildings_subtracted_3_after_CleanConditions".format(num_test)
building.WriteMdpaOutput(mdpa_out_name)

# we set the main_model_part to run CleanConditionsAngles
main_model_part = building.GetGeoModelPart()

stop_03 = time.time()
# we delete the condition in the angles
KratosGeo.CleaningUtilities(main_model_part).CleanConditionsAngles()
stop_04 = time.time()
print("[DEBUG PY] CleanConditionsAngles done in {} s".format(stop_04-stop_03))

# we set the geo model part to save the GiD file and mdpa file
building.SetGeoModelPart(main_model_part)

# we write the GiD file
building.CreateGidControlOutput("cfd_data/test_{}/gid_file/Ascii/11_Box_buildings_after_CleanConditionsAngles".format(num_test), "GiD_PostAscii")
building.CreateGidControlOutput("cfd_data/test_{}/gid_file/Binary/11_Box_buildings_after_CleanConditionsAngles".format(num_test), "GiD_PostBinary")
# writing file mdpa
mdpa_out_name = "cfd_data/test_{}/mdpa_file/11_Box_buildings_after_CleanConditionsAngles".format(num_test)
building.WriteMdpaOutput(mdpa_out_name)


#############
### C F D ###
#############
main_model_part = building.GetGeoModelPart()

stop_03 = time.time()
model = GeoModel()
model.SetGeoModelPart(main_model_part)
model.GenerateCfdModelPart()

stop_A = time.time()
KratosGeo.FillCfdModelpartUtilities(main_model_part).FillModelPart(model.GetGeoCfdModelPart())
print("[DEBUG PY]\t-> FillModelPart filled in ", (time.time()-stop_A))

# model.FillPartsFluid("Parts_Fluid")
# stop_03_1 = time.time()
# print("[DEBUG PY]\t-> Parts_Fluid filled in ", (stop_03_1-stop_03))

# model.FillNoslip("SKIN_ISOSURFACE")
# stop_03_2 = time.time()
# print("[DEBUG PY]\t-> SKIN_ISOSURFACE filled in ", (stop_03_2-stop_03_1))

# model.FillSlip("TopModelPart")
# stop_03_3 = time.time()
# print("[DEBUG PY]\t-> TopModelPart filled in ", (stop_03_3-stop_03_2))

# model.FillSlip("BottomModelPart")
# stop_03_4 = time.time()
# print("[DEBUG PY]\t-> BottomModelPart filled in ", (stop_03_4-stop_03_3))

# model.FillInlet("Inlet")
# stop_03_5 = time.time()
# print("[DEBUG PY]\t-> Inlet filled in ", (stop_03_5-stop_03_4))

# model.FillOutlet("Outlet")
# stop_04 = time.time()
# print("[DEBUG PY]\t-> Outlet filled in ", (stop_04-stop_03_5))

# print("\n[DEBUG PY] Filled time: {}\n".format(stop_04-stop_03))

print("\n[DEBUG PY] GeoCfdModelPart AFTER\n", model.GetGeoCfdModelPart())

# writing file mdpa
model.SetGeoModelPart(model.GetGeoCfdModelPart())
mdpa_out_name = "cfd_data/test_{}/analysis_file/12_CFD_model".format(num_test)
model.WriteMdpaOutput(mdpa_out_name)

print("[DEBUG PY] Time: ", time.time() - start_time)
print("\n[DEBUG PY] TEST ", num_test, "END\n\n")
