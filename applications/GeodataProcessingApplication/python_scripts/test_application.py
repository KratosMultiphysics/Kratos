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
num_test = "A_006"
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
extract_center = True

import_terrain_mdpa = False			# INPUT PARAMETER

if import_terrain_mdpa:
	print("\n*** IMPORT TERRAIN FROM MDPA FILE ***\n")

	# import domain from mdpa file
	importer._InitializeModelPart("test_model")
	terrain_model_part = importer.ModelPart
	model_part_in = "data/mdpa_file/domain_25_sept_2019/01_Mesh_cylinder"
	
	KratosMultiphysics.ModelPartIO(model_part_in).ReadModelPart(terrain_model_part)

	mesher.SetGeoModelPart(terrain_model_part)
	mesher.CreateGidControlOutput("cfd_data/test_{}/gid_file/Ascii/01_Mesh_cylinder".format(num_test), "GiD_PostAscii")
	mesher.CreateGidControlOutput("cfd_data/test_{}/gid_file/Binary/01_Mesh_cylinder".format(num_test), "GiD_PostBinary")

	# if there is CenterCondition we create center_model_part
	if (terrain_model_part.HasSubModelPart("CenterCondition")):
		# print("\n*** CenterCondition ***\n")
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
	model_part_in = "data/mdpa_file/domain_25_sept_2019/03_Mesh_cylinder_refinement_1"

	KratosMultiphysics.ModelPartIO(model_part_in).ReadModelPart(terrain_model_part)

	mesher.SetGeoModelPart(terrain_model_part)
	mesher.CreateGidControlOutput("cfd_data/test_{}/gid_file/Ascii/03_Mesh_cylinder_refinement_1".format(num_test), "GiD_PostAscii")
	mesher.CreateGidControlOutput("cfd_data/test_{}/gid_file/Binary/03_Mesh_cylinder_refinement_1".format(num_test), "GiD_PostBinary")

else:
	print("\n*** IMPORT TERRAIN FROM STL FILE ***\n")

	# import STL terrain and compute mesh circle
	stl_name = "terrain_2"
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
	height = 100.0		# height = 80.0
	# elem_list = mesher.MeshCircleWithTerrainPoints(height, 40, True)	# if the value is True we can save a list with the elements that are inside r_buildings
	
	# print("MeshCircleWithTerrainPoints")
	# mesher.MeshCircleWithTerrainPoints(height, 40, extract_center)	# if extract_center = True we create a sub model part with conditions that are inside r_buildings
	
	num_sector = 37
	print("MeshCircleWithTerrainPoints_old")
	mesher.MeshCircleWithTerrainPoints_old(height, 40, num_sector, extract_center)	# if extract_center = True we create a sub model part with conditions that are inside r_buildings

	# we populate the sub model part sectors
	mesher.MeshSectors(num_sector)
	# for n_sect in range(1, 13):
	# 	current_smp = mesher.ModelPart.GetSubModelPart("LateralSector_{}".format(n_sect))
	# 	print("LateralSector_{}".format(n_sect))
	# 	for cond in current_smp.Conditions:
	# 		print("\t", cond.Id)
	# 		for node in cond.GetNodes():
	# 			print("\t\t", node.Id)
	
	### creo dei model part nuovi per poterli visualizzare in GiD ###
	for n_sect in range(1, num_sector+1):
		# creo il nuovo model part che devo visualizzare
		importer._InitializeModelPart("LateralSector_{}".format(n_sect))
		new_sect_mp = importer.ModelPart
		properties = new_sect_mp.Properties[0]
		# leggo il smp del settore n-esimo
		current_smp = mesher.ModelPart.GetSubModelPart("LateralSector_{}".format(n_sect))

		# creo i nodi nel model part (copiando le informazioni dal smp)
		for cond in current_smp.Conditions:
			nodes = []
			for node in cond.GetNodes():
				nodes.append(node.Id)
				new_sect_mp.CreateNewNode(node.Id, node.X, node.Y, node.Z)
			new_sect_mp.CreateNewCondition("SurfaceCondition3D3N", cond.Id, nodes, properties)
		
		mesher_sectors = GeoMesher()
		mesher_sectors.SetGeoModelPart(new_sect_mp)
		mesher_sectors.CreateGidControlOutput("cfd_data/test_{}/gid_file/Ascii/LateralSector_{}".format(num_test, n_sect), "GiD_PostAscii")
		mesher_sectors.CreateGidControlOutput("cfd_data/test_{}/gid_file/Binary/LateralSector_{}".format(num_test, n_sect), "GiD_PostBinary")
	### FINE ###

	# writing GiD file
	mesher.CreateGidControlOutput("cfd_data/test_{}/gid_file/Ascii/01_Mesh_cylinder".format(num_test), "GiD_PostAscii")
	mesher.CreateGidControlOutput("cfd_data/test_{}/gid_file/Binary/01_Mesh_cylinder".format(num_test), "GiD_PostBinary")
	input("STOOOOOP")

	###################################################################
	# we split "LateralModelPart" into "Inlet" and "Outlet". We remove "LateralModelPart" after
	model_part = mesher.GetGeoModelPart()
	inlet_nodes = [];	inlet_cond = []
	outlet_nodes =[];	outlet_cond = []
	for cond in model_part.GetSubModelPart("LateralModelPart").Conditions:
		nodes = cond.GetNodes()
		if (nodes[0].X <= 0) and (nodes[1].X <= 0) and (nodes[2].X <= 0):
			inlet_nodes.extend([nodes[0].Id, nodes[1].Id, nodes[2].Id])
			inlet_cond.append(cond.Id)
		else:
			outlet_nodes.extend([nodes[0].Id, nodes[1].Id, nodes[2].Id])
			outlet_cond.append(cond.Id)
	
	inlet_model_part = model_part.CreateSubModelPart("Inlet")
	inlet_model_part.AddNodes(inlet_nodes)
	inlet_model_part.AddConditions(inlet_cond)

	outlet_model_part = model_part.CreateSubModelPart("Outlet")
	outlet_model_part.AddNodes(outlet_nodes)
	outlet_model_part.AddConditions(outlet_cond)

	model_part.RemoveSubModelPart("LateralModelPart")
	###################################################################

	terrain_model_part = mesher.GetGeoModelPart()
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

	# input("STOP 1")

	### MMG ###
	# 1st ground refinement
	print("\n\n***** 1st ground refinement *****\n")
	# distance field from ground
	mesher.ComputeDistanceFieldFromGround()
	mesher.CreateGidControlOutput("cfd_data/test_{}/gid_file/Ascii/02_Mesh_cylinder_distance_field_1".format(num_test), "GiD_PostAscii")
	mesher.CreateGidControlOutput("cfd_data/test_{}/gid_file/Binary/02_Mesh_cylinder_distance_field_1".format(num_test), "GiD_PostBinary")
	# writing file mdpa
	mdpa_out_name = "cfd_data/test_{}/mdpa_file/02_Mesh_cylinder_distance_field_1".format(num_test)
	mesher.WriteMdpaOutput(mdpa_out_name)
	# refine mesh
	# mesher.RefineMesh_test(5.0, 10.0)
	mesher.RefineMesh_test(10.0, 100.0)
	# mesher.RefineMesh_test(50.0, 200.0)
	mesher.CreateGidControlOutput("cfd_data/test_{}/gid_file/Ascii/03_Mesh_cylinder_refinement_1".format(num_test), "GiD_PostAscii")
	mesher.CreateGidControlOutput("cfd_data/test_{}/gid_file/Binary/03_Mesh_cylinder_refinement_1".format(num_test), "GiD_PostBinary")
	# writing file mdpa
	mdpa_out_name = "cfd_data/test_{}/mdpa_file/03_Mesh_cylinder_refinement_1".format(num_test)
	mesher.WriteMdpaOutput(mdpa_out_name)
	print("BOX REFINEMENT 1 DONE!")
	# input("STOP 2")

### creo dei model part nuovi per poterli visualizzare in GiD ###
for n_sect in range(1, 13):
	# creo il nuovo model part che devo visualizzare
	importer._InitializeModelPart("LateralSector_{}".format(n_sect))
	new_sect_mp = importer.ModelPart
	properties = new_sect_mp.Properties[0]
	# leggo il smp del settore n-esimo
	current_smp = mesher.ModelPart.GetSubModelPart("LateralSector_{}".format(n_sect))

	# creo i nodi nel model part (copiando le informazioni dal smp)
	for cond in current_smp.Conditions:
		nodes = []
		for node in cond.GetNodes():
			nodes.append(node.Id)
			new_sect_mp.CreateNewNode(node.Id, node.X, node.Y, node.Z)
		new_sect_mp.CreateNewCondition("SurfaceCondition3D3N", cond.Id, nodes, properties)
	
	mesher_sectors = GeoMesher()
	mesher_sectors.SetGeoModelPart(new_sect_mp)
	mesher_sectors.CreateGidControlOutput("cfd_data/test_{}/gid_file/Ascii/LateralSector_{}".format(num_test, n_sect), "GiD_PostAscii")
	mesher_sectors.CreateGidControlOutput("cfd_data/test_{}/gid_file/Binary/LateralSector_{}".format(num_test, n_sect), "GiD_PostBinary")
### FINE ###

print("\n** TERRAIN DONE! **\n")
input("PAUSE (press to continue)")

# terrain volume mesh ready
main_model_part = mesher.GetGeoModelPart()


# only for check (delete it)
print(center_model_part)
mesher.SetGeoModelPart(center_model_part)
mesher.CreateGidControlOutput("cfd_data/test_{}/gid_file/Ascii/center_model_part".format(num_test), "GiD_PostAscii")
mesher.CreateGidControlOutput("cfd_data/test_{}/gid_file/Binary/center_model_part".format(num_test), "GiD_PostBinary")

building = GeoBuilding()
building.SetGeoModelPart(main_model_part)

# import buildings
import_building_mdpa = False
if import_building_mdpa:
	print("\n*** IMPORT BUILDINGS FROM MDPA FILE ***\n")
	# import domain from mdpa file
	model_part_in = "/data/mdpa_file/buildings_num_20"
	building.ImportBuildingHullMDPA(model_part_in)
	building_model_part = building.GetBuildingModelPart()

	# importer._InitializeModelPart("BuildingModelPart")
	# building_model_part = importer.ModelPart
	# model_part_in = "data/mdpa_file/buildings_num_9"	# model_part_in = "data/mdpa_file/buildings_100"
	# KratosMultiphysics.ModelPartIO(model_part_in).ReadModelPart(building_model_part)
	# # importer.SetGeoModelPart(building_model_part)
else:
	print("\n*** IMPORT BUILDINGS FROM OBJ FILE ***\n")
	# obj file name
	obj_file_in = "data/buildings/buildings.obj"
	importer.ObjImportBuildings(obj_file_in, "BuildingModelPart", True)		# (obj_file_name_input, name_model_part, change_coord)
	building_model_part = importer.GetGeoModelPart()

# we write GiD file
mesher.SetGeoModelPart(building_model_part)
mesher.CreateGidControlOutput("cfd_data/test_{}/gid_file/Ascii/buildings".format(num_test), "GiD_PostAscii")
mesher.CreateGidControlOutput("cfd_data/test_{}/gid_file/Binary/buildings".format(num_test), "GiD_PostBinary")

# writing file mdpa
mdpa_out_name = "cfd_data/test_{}/mdpa_file/buildings".format(num_test)
mesher.WriteMdpaOutput(mdpa_out_name)

input("\nIMPORT BUILDINGS DONE! (press to continue)")

# # FOR DEBUG PURPOSE ONLY
# STL_center_name = "cfd_data/test_{}/stl_file/terrain_center.stl".format(num_test)
# preprocessor.WriteSTL(STL_center_name, elem_list)
# importer.StlImport(STL_center_name)

if extract_center:
	building.SetGeoModelPart(building_model_part)
	building.ShiftBuildingOnTerrain( center_model_part)
	# building.ShiftBuildingOnTerrain(building_model_part, center_model_part)

	# delete building placed out of center_model_part. If a building is out of range of the center_model_part is placed on z = 0
	# building.DeleteBuildingsUnderValue(building_model_part, 10)		# TODO: CHECK IT
	building.CreateGidControlOutput("cfd_data/test_{}/gid_file/Ascii/buildings_shifted".format(num_test), "GiD_PostAscii")
	building.CreateGidControlOutput("cfd_data/test_{}/gid_file/Binary/buildings_shifted".format(num_test), "GiD_PostBinary")
	# writing file mdpa
	mdpa_out_name = "cfd_data/test_{}/mdpa_file/buildings_shifted".format(num_test)
	building.WriteMdpaOutput(mdpa_out_name)
input("\nPAUSE shift buildings (press to continue)")


##############################################################################
# building steps
print("\n\nSTART PROCESS BUILDING")
stop_begin_loop = time.time()

# building = GeoBuilding()
building.SetGeoModelPart(main_model_part)		# it is a duplicate. CHECK IT
building.ImportBuilding(building_model_part)

# # 1st cut
# print("\n\n***** STEP 1 *****\n")
# # distance field from hull
# building.ComputeDistanceFieldFromHull(False, 1e-7)
# building.CreateGidControlOutput("cfd_data/test_{}/gid_file/Ascii/04_Box_buildings_distance_field_1".format(num_test), "GiD_PostAscii")
# building.CreateGidControlOutput("cfd_data/test_{}/gid_file/Binary/04_Box_buildings_distance_field_1".format(num_test), "GiD_PostBinary")
# # writing file mdpa
# mdpa_out_name = "cfd_data/test_{}/mdpa_file/04_Box_buildings_distance_field_1".format(num_test)
# building.WriteMdpaOutput(mdpa_out_name)
# # input("\n04_Box_buildings_distance_field_1 DONE! (press to continue)")
# # subtract buildings
# # building.SubtractBuildingMOD(0.5, 10.0, 0.1, "Linear", "STANDARD", "false")			# interpolation = constant; disc_type = STANDARD; remove_internal_regions=false
# building.SubtractBuildingMOD(5.0, 100.0, 0.1, "Linear", "STANDARD", "false")			# interpolation = constant; disc_type = STANDARD; remove_internal_regions=false
# building.CreateGidControlOutput("cfd_data/test_{}/gid_file/Ascii/05_Box_buildings_subtracted_1".format(num_test), "GiD_PostAscii")
# building.CreateGidControlOutput("cfd_data/test_{}/gid_file/Binary/05_Box_buildings_subtracted_1".format(num_test), "GiD_PostBinary")
# # writing file mdpa
# mdpa_out_name = "cfd_data/test_{}/mdpa_file/05_Box_buildings_subtracted_1".format(num_test)
# building.WriteMdpaOutput(mdpa_out_name)
# # input("\n05_Box_buildings_subtracted_1 DONE! (press to continue)")

# # 2nd cut
# print("\n\n***** STEP 2 *****\n")
# # distance field from hull
# building.ComputeDistanceFieldFromHull(False, 1e-7)
# building.CreateGidControlOutput("cfd_data/test_{}/gid_file/Ascii/06_Box_buildings_distance_field_2".format(num_test), "GiD_PostAscii")
# building.CreateGidControlOutput("cfd_data/test_{}/gid_file/Binary/06_Box_buildings_distance_field_2".format(num_test), "GiD_PostBinary")
# # writing file mdpa
# mdpa_out_name = "cfd_data/test_{}/mdpa_file/06_Box_buildings_distance_field_2".format(num_test)
# building.WriteMdpaOutput(mdpa_out_name)
# # subtract buildings
# # building.SubtractBuildingMOD(0.2, 10.0, 0.1, "constant", "STANDARD", "false")			# interpolation = linear; disc_type = STANDARD; remove_internal_regions=false
# building.SubtractBuildingMOD(3.0, 10.0, 0.1, "Constant", "STANDARD", "false")			# interpolation = linear; disc_type = STANDARD; remove_internal_regions=false
# building.CreateGidControlOutput("cfd_data/test_{}/gid_file/Ascii/07_Box_buildings_subtracted_2".format(num_test), "GiD_PostAscii")
# building.CreateGidControlOutput("cfd_data/test_{}/gid_file/Binary/07_Box_buildings_subtracted_2".format(num_test), "GiD_PostBinary")
# # writing file mdpa
# mdpa_out_name = "cfd_data/test_{}/mdpa_file/07_Box_buildings_subtracted_2".format(num_test)
# building.WriteMdpaOutput(mdpa_out_name)

# 3rd cut
print("\n\n***** STEP 3 *****\n")
# distance field from hull
building.ComputeDistanceFieldFromHull(False, 1e-7)
building.CreateGidControlOutput("cfd_data/test_{}/gid_file/Ascii/08_Box_buildings_distance_field_3".format(num_test), "GiD_PostAscii")
building.CreateGidControlOutput("cfd_data/test_{}/gid_file/Binary/08_Box_buildings_distance_field_3".format(num_test), "GiD_PostBinary")
# writing file mdpa
mdpa_out_name = "cfd_data/test_{}/mdpa_file/08_Box_buildings_distance_field_3".format(num_test)
building.WriteMdpaOutput(mdpa_out_name)
# input("\n08_Box_buildings_distance_field_3 (press to continue)")
# subtract buildings
# building.SubtractBuildingMOD(0.2, 10.0, 0.1, "exponential", "ISOSURFACE", "true")			# interpolation = exponential; disc_type = ISOSURFACE; remove_internal_regions=true
building.SubtractBuildingMOD(2.0, 100.0, 0.1, "exponential", "ISOSURFACE", "true")			# interpolation = exponential; disc_type = ISOSURFACE; remove_internal_regions=true
building.CreateGidControlOutput("cfd_data/test_{}/gid_file/Ascii/09_Box_buildings_subtracted_3_before_CleanConditions".format(num_test), "GiD_PostAscii")
building.CreateGidControlOutput("cfd_data/test_{}/gid_file/Binary/09_Box_buildings_subtracted_3_before_CleanConditions".format(num_test), "GiD_PostBinary")

# we update main_model_part
main_model_part = building.GetGeoModelPart()

stop_01 = time.time()
# we delete the condition if at least one node it is not in main model part
KratosGeo.CleaningUtilities(main_model_part).CleanConditions()
stop_02 = time.time()
print("CleanConditions done in ", (stop_02 - stop_01))

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
print("CleanConditionsAngles done in {} s".format(stop_04-stop_03))

# we set the geo model part to save the GiD file and mdpa file
building.SetGeoModelPart(main_model_part)

# we write the GiD file
building.CreateGidControlOutput("cfd_data/test_{}/gid_file/Ascii/11_Box_buildings_after_CleanConditionsAngles".format(num_test), "GiD_PostAscii")
building.CreateGidControlOutput("cfd_data/test_{}/gid_file/Binary/11_Box_buildings_after_CleanConditionsAngles".format(num_test), "GiD_PostBinary")
# writing file mdpa
mdpa_out_name = "cfd_data/test_{}/mdpa_file/11_Box_buildings_after_CleanConditionsAngles".format(num_test)
building.WriteMdpaOutput(mdpa_out_name)





##################################################
# CFD PART
main_model_part = building.GetGeoModelPart()

stop_03 = time.time()
model = GeoModel()
model.SetGeoModelPart(main_model_part)
model.GenerateCfdModelPart()

model.FillPartsFluid("Parts_Fluid")
stop_03_1 = time.time()
print("\t-> Parts_Fluid filled in ", (stop_03_1-stop_03))

model.FillNoslip("SKIN_ISOSURFACE")
stop_03_2 = time.time()
print("\t-> SKIN_ISOSURFACE filled in ", (stop_03_2-stop_03_1))

model.FillSlip("TopModelPart")
stop_03_3 = time.time()
print("\t-> TopModelPart filled in ", (stop_03_3-stop_03_2))

model.FillSlip("BottomModelPart")
stop_03_4 = time.time()
print("\t-> BottomModelPart filled in ", (stop_03_4-stop_03_3))

model.FillInlet("Inlet")
stop_03_5 = time.time()
print("\t-> Inlet filled in ", (stop_03_5-stop_03_4))

model.FillOutlet("Outlet")
stop_04 = time.time()
print("\t-> Inlet filled in ", (stop_04-stop_03_5))

print("\n*** Filled time: {} ***\n".format(stop_04-stop_03))

print("\nGeoCfdModelPart AFTER\n", model.GetGeoCfdModelPart())

# writing file mdpa
model.SetGeoModelPart(model.GetGeoCfdModelPart())
mdpa_out_name = "cfd_data/test_{}/analysis_file/12_CFD_model".format(num_test)
model.WriteMdpaOutput(mdpa_out_name)


##################################################
# JSON file
print("\n\n*** JSON part ***")
stop_05 = time.time()
# we copy MainKratos.py file
from shutil import copyfile
src_path = "data/MainKratos.py"
dest_path = "cfd_data/test_{}/analysis_file/MainKratos.py".format(num_test)
copyfile(src_path, dest_path)

# names of the all parts
inlet_name = "Inlet"
outlet_name = "Outlet"
slip_name = ["TopModelPart", "BottomModelPart"]		# it can be a string or a list of strings
noslip_name = "SKIN_ISOSURFACE"				# it can be a string or a list of strings
volume_name = "Parts_Fluid"

# we read the JSON file (the base form) for the CFD analysis
input_string = open("data/parameters/test/ProjectParameters.json",'r').read()
settings = KratosMultiphysics.Parameters(input_string)

# we set "problem_name"
problem_name = "test_{}".format(num_test)
settings["problem_data"].AddEmptyValue("problem_name")
settings["problem_data"]["problem_name"].SetString(problem_name)

# we set "output_name"
gid_output = settings["output_processes"]["gid_output"]		# gid_output is an array with only one element (where there all the entries are)
gid_output[0]["Parameters"].AddEmptyValue("output_name")
gid_output[0]["Parameters"]["output_name"].SetString(problem_name)

# we set "model_import_settings"
settings["solver_settings"]["model_import_settings"].AddEmptyValue("input_filename")
settings["solver_settings"]["model_import_settings"]["input_filename"].SetString(mdpa_out_name)		# the name of the last model part saved

# we set "volume_model_part_name"
settings["solver_settings"].AddEmptyValue("volume_model_part_name")
settings["solver_settings"]["volume_model_part_name"].SetString(volume_name)

# we set "skin_parts"
if (isinstance(slip_name, list)):
	settings["solver_settings"].RemoveValue("skin_parts")	# to avoid duplicates
	settings["solver_settings"].AddEmptyArray("skin_parts")
	for i in range(len(slip_name)):
		settings["solver_settings"]["skin_parts"].Append(slip_name[i])

elif (isinstance(slip_name, str)):
	settings["solver_settings"]["skin_parts"].SetString(slip_name)

# we set "no_skin_parts" after a type check
if (isinstance(noslip_name, list)):
	settings["solver_settings"].RemoveValue("no_skin_parts")	# to avoid duplicates
	settings["solver_settings"].AddEmptyArray("no_skin_parts")
	for i in range(len(noslip_name)):
		settings["solver_settings"]["no_skin_parts"].Append(noslip_name[i])

elif (isinstance(noslip_name, str)):
	settings["solver_settings"]["no_skin_parts"].SetString(noslip_name)

# we set "boundary_conditions_process_list" array
settings["processes"].RemoveValue("boundary_conditions_process_list")	# to avoid duplicates
settings["processes"].AddEmptyArray("boundary_conditions_process_list")

# inlet
inlet_process_1 = KratosMultiphysics.Parameters("""
				{
					"python_module" : "apply_inlet_process",
					"kratos_module" : "KratosMultiphysics.FluidDynamicsApplication",
					"Parameters"    : {
						"model_part_name" : "FluidModelPart.""" + str(inlet_name) + """",
						"variable_name"   : "VELOCITY",
						"modulus"         : 1.0,
						"direction"       : "automatic_inwards_normal",
						"interval"        : [0,1]
					}
				} """)
settings["processes"]["boundary_conditions_process_list"].Append(inlet_process_1)

inlet_process_2 = KratosMultiphysics.Parameters("""
				{
					"python_module" : "apply_inlet_process",
					"kratos_module" : "KratosMultiphysics.FluidDynamicsApplication",
					"Parameters"    : {
						"model_part_name" : "FluidModelPart.""" + str(inlet_name) + """",
						"variable_name"   : "VELOCITY",
						"modulus"         : 0.5,
						"direction"       : "automatic_inwards_normal",
						"interval"        : [1,"End"]
					}
				} """)
settings["processes"]["boundary_conditions_process_list"].Append(inlet_process_2)

# outlet
outlet_process = KratosMultiphysics.Parameters("""
				{
					"python_module" : "apply_outlet_process",
					"kratos_module" : "KratosMultiphysics.FluidDynamicsApplication",
					"Parameters"    : {
						"model_part_name" : "FluidModelPart.""" + str(outlet_name) + """",
						"variable_name"   : "PRESSURE",
						"constrained"        : true,
						"value"              : 0.0,
						"hydrostatic_outlet" : false,
						"h_top"              : 0.0
					}
				} """)
settings["processes"]["boundary_conditions_process_list"].Append(outlet_process)

# slip (it can be a string or a list of strings)
if (isinstance(slip_name, list)):
	for i_slip in slip_name:
		slip_process = KratosMultiphysics.Parameters("""
						{
							"python_module" : "apply_slip_process",
							"kratos_module" : "KratosMultiphysics.FluidDynamicsApplication",
							"Parameters"    : {
								"model_part_name" : "FluidModelPart.""" + str(i_slip) + """"
							}
						} """)
		settings["processes"]["boundary_conditions_process_list"].Append(slip_process)
elif (isinstance(slip_name, str)):
	slip_process = KratosMultiphysics.Parameters("""
					{
						"python_module" : "apply_slip_process",
						"kratos_module" : "KratosMultiphysics.FluidDynamicsApplication",
						"Parameters"    : {
							"model_part_name" : "FluidModelPart.""" + str(slip_name) + """"
						}
					} """)
	settings["processes"]["boundary_conditions_process_list"].Append(slip_process)

# noslip (it can be a string or a list of strings)
if (isinstance(noslip_name, list)):
	for i_noslip in noslip_name:
		noslip_process = KratosMultiphysics.Parameters("""
						{
							"python_module" : "apply_noslip_process",
							"kratos_module" : "KratosMultiphysics.FluidDynamicsApplication",
							"Parameters"    : {
								"model_part_name" : "FluidModelPart.""" + str(i_noslip) + """"
							}
						} """)
		settings["processes"]["boundary_conditions_process_list"].Append(noslip_process)
elif (isinstance(noslip_name, str)):
	noslip_process = KratosMultiphysics.Parameters("""
					{
						"python_module" : "apply_noslip_process",
						"kratos_module" : "KratosMultiphysics.FluidDynamicsApplication",
						"Parameters"    : {
							"model_part_name" : "FluidModelPart.""" + str(noslip_name) + """"
						}
					} """)
	settings["processes"]["boundary_conditions_process_list"].Append(noslip_process)

# we save the modified json file
file_name = "cfd_data/test_{}/analysis_file/ProjectParameters.json".format(num_test)
with open(file_name, mode="w", encoding="utf-8") as file_json:
	file_json.write(settings.PrettyPrintJsonString())

stop_06 = time.time()
print("\n*** JSON time: {} ***\n".format(stop_06-stop_05))
# end JSON file
##################################################

print("*** Time: ", time.time() - start_time)
print("\nTEST ", num_test, "END\n\n")
