""" cylindrical domain
	#
	# create:	29 August 2019 (copy from test_cylinder_14_import_mdpa.py)
	# edit:		30 August 2019
	# edit:		02 September 2019
	# edit:		03 September 2019 -> added CreateGidControlOutput with GiD_PostAscii and GiD_PostBinary
	# edit:		04 September 2019 -> building MDPA file are imported with ImportBuildingHullMDPA function now
	# edit:		05 September 2019 -> in import_terrain_mdpa=False the element now is "Element2D3N" (in center_model_part)
	# edit:		09 September 2019 -> added process to fill Parts_Fluid, Inlet, Outlet and Slip sub model part

	# we can import the domain mdpa file and building mdpa;
	# otherwise we can import terrain in STL file and building in OBJ file respectively (and run MMG process).
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
num_test = "28_september_20buildings_2_refinements"
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

import_terrain_mdpa = True
if import_terrain_mdpa:
	print("\n*** IMPORT TERRAIN FROM MDPA FILE ***\n")

	# import domain from mdpa file
	importer._InitializeModelPart("test_model")
	terrain_model_part = importer.ModelPart
	# model_part_in = "data/mdpa_file/01_Mesh_cylinder"
	# model_part_in = "data/mdpa_file/domain_06_sept_2019/01_Mesh_cylinder"
	model_part_in = "data/mdpa_file/domain_09_sept_2019/01_Mesh_cylinder"
	
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
			# center_model_part.CreateNewCondition("SurfaceCondition3D3N", cond.Id, [nodes[0].Id, nodes[1].Id, nodes[2].Id], prop)
			# center_model_part.CreateNewElement("Element3D3N", cond.Id, [nodes[0].Id, nodes[1].Id, nodes[2].Id], prop)
			center_model_part.CreateNewElement("Element2D3N", cond.Id, [nodes[0].Id, nodes[1].Id, nodes[2].Id], prop)	# we need this to shift the buildings on terrain
	
	# we clear the model part to import the domain after the MMG process
	importer._InitializeModelPart("test_model_after_MMG")
	terrain_model_part = importer.ModelPart
	# model_part_in = "data/mdpa_file/03_Mesh_cylinder_refinement_1"
	# model_part_in = "data/mdpa_file/domain_06_sept_2019/03_Mesh_cylinder_refinement_1"
	model_part_in = "data/mdpa_file/domain_09_sept_2019/03_Mesh_cylinder_refinement_1"

	KratosMultiphysics.ModelPartIO(model_part_in).ReadModelPart(terrain_model_part)

	mesher.SetGeoModelPart(terrain_model_part)
	mesher.CreateGidControlOutput("cfd_data/test_{}/gid_file/Ascii/03_Mesh_cylinder_refinement_1".format(num_test), "GiD_PostAscii")
	mesher.CreateGidControlOutput("cfd_data/test_{}/gid_file/Binary/03_Mesh_cylinder_refinement_1".format(num_test), "GiD_PostBinary")

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
	
	# print("MeshCircleWithTerrainPoints")
	# mesher.MeshCircleWithTerrainPoints(height, 40, extract_center)	# if extract_center = True we create a sub model part with conditions that are inside r_buildings
	
	print("MeshCircleWithTerrainPoints_old")
	mesher.MeshCircleWithTerrainPoints_old(height, 40, extract_center)	# if extract_center = True we create a sub model part with conditions that are inside r_buildings

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
			# center_model_part.CreateNewCondition("SurfaceCondition3D3N", cond.Id, [nodes[0].Id, nodes[1].Id, nodes[2].Id], prop)
			# center_model_part.CreateNewElement("Element3D3N", cond.Id, [nodes[0].Id, nodes[1].Id, nodes[2].Id], prop)
			center_model_part.CreateNewElement("Element2D3N", cond.Id, [nodes[0].Id, nodes[1].Id, nodes[2].Id], prop)
		
	# we need to remove this sub model part to avoid problems in MMG process
	mesher.GetGeoModelPart().RemoveSubModelPart("CenterCondition")

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
	# mesher.RefineMesh_test(10.0, 100.0)
	mesher.RefineMesh_test(50.0, 200.0)
	mesher.CreateGidControlOutput("cfd_data/test_{}/gid_file/Ascii/03_Mesh_cylinder_refinement_1".format(num_test), "GiD_PostAscii")
	mesher.CreateGidControlOutput("cfd_data/test_{}/gid_file/Binary/03_Mesh_cylinder_refinement_1".format(num_test), "GiD_PostBinary")
	# writing file mdpa
	mdpa_out_name = "cfd_data/test_{}/mdpa_file/03_Mesh_cylinder_refinement_1".format(num_test)
	mesher.WriteMdpaOutput(mdpa_out_name)
	print("BOX REFINEMENT 1 DONE!")
	# input("PAUSE")

print("\n** TERRAIN DONE! **\n")
# input("PAUSE (press to continue)")

# terrain volume mesh ready
main_model_part = mesher.GetGeoModelPart()


# only for check (delete it)
print(center_model_part)
mesher.SetGeoModelPart(center_model_part)
mesher.CreateGidControlOutput("cfd_data/test_{}/gid_file/Ascii/center_model_part".format(num_test), "GiD_PostAscii")
mesher.CreateGidControlOutput("cfd_data/test_{}/gid_file/Binary/center_model_part".format(num_test), "GiD_PostBinary")
# input("\n03_Mesh_cylinder_refinement_1 DONE!\n")

building = GeoBuilding()
building.SetGeoModelPart(main_model_part)

# import buildings
import_building_mdpa = True
if import_building_mdpa:
	print("\n*** IMPORT BUILDINGS FROM MDPA FILE ***\n")
	# import domain from mdpa file
	# model_part_in = "/data/mdpa_file/buildings_num_9"	# model_part_in = "data/mdpa_file/buildings_100"
	# model_part_in = "/data/mdpa_file/buildings_num_6"	# model_part_in = "data/mdpa_file/buildings_100"
	model_part_in = "/data/mdpa_file/buildings_num_20"	# model_part_in = "data/mdpa_file/buildings_100"
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
	obj_file_in = "data/buildings/buildings_barcellona_1.obj"
	importer.ObjImport(obj_file_in, "BuildingModelPart", False)		# with "False" we avoid the change of the coordinates
	building_model_part = importer.GetGeoModelPart()

# we write GiD file
mesher.SetGeoModelPart(building_model_part)
mesher.CreateGidControlOutput("cfd_data/test_{}/gid_file/Ascii/04_buildings".format(num_test), "GiD_PostAscii")
mesher.CreateGidControlOutput("cfd_data/test_{}/gid_file/Binary/04_buildings".format(num_test), "GiD_PostBinary")

# input("\nIMPORT BUILDINGS DONE! (press to continue)")

# # FOR DEBUG PURPOSE ONLY
# STL_center_name = "cfd_data/test_{}/stl_file/terrain_center.stl".format(num_test)
# preprocessor.WriteSTL(STL_center_name, elem_list)
# importer.StlImport(STL_center_name)


if extract_center:
	building.ShiftBuildingOnTerrain(building_model_part, center_model_part)
	# delete building placed out of center_model_part. If a building is out of range of the center_model_part is placed on z = 0
	building.DeleteBuildingsUnderValue(building_model_part, 10)		# CHECK IT
	building.CreateGidControlOutput("cfd_data/test_{}/gid_file/Ascii/buildings_shifted".format(num_test), "GiD_PostAscii")
	building.CreateGidControlOutput("cfd_data/test_{}/gid_file/Binary/buildings_shifted".format(num_test), "GiD_PostBinary")
	# writing file mdpa
	mdpa_out_name = "cfd_data/test_{}/mdpa_file/buildings_shifted".format(num_test)
	building.WriteMdpaOutput(mdpa_out_name)
# input("\nPAUSE shift buildings (press to continue)")


##############################################################################
# building steps
print("\n\nSTART PROCESS BUILDING")
stop_begin_loop = time.time()

# building = GeoBuilding()
building.SetGeoModelPart(main_model_part)		# it is a duplicate. CHECK IT
building.ImportBuilding(building_model_part)

# 1st cut
print("\n\n***** STEP 1 *****\n")
# distance field from hull
building.ComputeDistanceFieldFromHull(False, 1e-7)
building.CreateGidControlOutput("cfd_data/test_{}/gid_file/Ascii/04_Box_buildings_distance_field_1".format(num_test), "GiD_PostAscii")
building.CreateGidControlOutput("cfd_data/test_{}/gid_file/Binary/04_Box_buildings_distance_field_1".format(num_test), "GiD_PostBinary")
# writing file mdpa
mdpa_out_name = "cfd_data/test_{}/mdpa_file/04_Box_buildings_distance_field_1".format(num_test)
building.WriteMdpaOutput(mdpa_out_name)
# input("\n04_Box_buildings_distance_field_1 DONE! (press to continue)")
# subtract buildings
# building.SubtractBuildingMOD(0.5, 10.0, 0.1, "Linear", "STANDARD", "false")			# interpolation = constant; disc_type = STANDARD; remove_internal_regions=false
building.SubtractBuildingMOD(5.0, 100.0, 0.1, "Linear", "STANDARD", "false")			# interpolation = constant; disc_type = STANDARD; remove_internal_regions=false
building.CreateGidControlOutput("cfd_data/test_{}/gid_file/Ascii/05_Box_buildings_subtracted_1".format(num_test), "GiD_PostAscii")
building.CreateGidControlOutput("cfd_data/test_{}/gid_file/Binary/05_Box_buildings_subtracted_1".format(num_test), "GiD_PostBinary")
# writing file mdpa
mdpa_out_name = "cfd_data/test_{}/mdpa_file/05_Box_buildings_subtracted_1".format(num_test)
building.WriteMdpaOutput(mdpa_out_name)
# input("\n05_Box_buildings_subtracted_1 DONE! (press to continue)")

# 2nd cut
print("\n\n***** STEP 2 *****\n")
# distance field from hull
building.ComputeDistanceFieldFromHull(False, 1e-7)
building.CreateGidControlOutput("cfd_data/test_{}/gid_file/Ascii/06_Box_buildings_distance_field_2".format(num_test), "GiD_PostAscii")
building.CreateGidControlOutput("cfd_data/test_{}/gid_file/Binary/06_Box_buildings_distance_field_2".format(num_test), "GiD_PostBinary")
# writing file mdpa
mdpa_out_name = "cfd_data/test_{}/mdpa_file/06_Box_buildings_distance_field_2".format(num_test)
building.WriteMdpaOutput(mdpa_out_name)
# subtract buildings
# building.SubtractBuildingMOD(0.2, 10.0, 0.1, "constant", "STANDARD", "false")			# interpolation = linear; disc_type = STANDARD; remove_internal_regions=false
building.SubtractBuildingMOD(3.0, 10.0, 0.1, "Constant", "STANDARD", "false")			# interpolation = linear; disc_type = STANDARD; remove_internal_regions=false
building.CreateGidControlOutput("cfd_data/test_{}/gid_file/Ascii/07_Box_buildings_subtracted_2".format(num_test), "GiD_PostAscii")
building.CreateGidControlOutput("cfd_data/test_{}/gid_file/Binary/07_Box_buildings_subtracted_2".format(num_test), "GiD_PostBinary")
# writing file mdpa
mdpa_out_name = "cfd_data/test_{}/mdpa_file/07_Box_buildings_subtracted_2".format(num_test)
building.WriteMdpaOutput(mdpa_out_name)

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
# writing file mdpa
mdpa_out_name = "cfd_data/test_{}/mdpa_file/09_Box_buildings_subtracted_3_before_CleanConditions".format(num_test)
building.WriteMdpaOutput(mdpa_out_name)
# input("\n09_Box_buildings_subtracted_3 (press to continue)")

# we update main_model_part
main_model_part = building.GetGeoModelPart()

# we delete the condition if at least one node it is not in main model part
KratosGeo.CleaningUtilities(main_model_part).CleanConditions()

# # we fill the bottom model part with the new conditions after MMG process
# KratosGeo.CleaningUtilities(main_model_part).FillBottom()

# we write gid file
building.CreateGidControlOutput("cfd_data/test_{}/gid_file/Ascii/10_Box_buildings_subtracted_3_after_CleanConditions".format(num_test), "GiD_PostAscii")
building.CreateGidControlOutput("cfd_data/test_{}/gid_file/Binary/10_Box_buildings_subtracted_3_after_CleanConditions".format(num_test), "GiD_PostBinary")
# we write mdpa file
mdpa_out_name = "cfd_data/test_{}/mdpa_file/10_Box_buildings_subtracted_3_after_CleanConditions".format(num_test)
building.WriteMdpaOutput(mdpa_out_name)



##################################################
# TODO: to add this part in geo_model.py

# CFD PART
main_model_part = building.GetGeoModelPart()

# we set the DENSITY and DYNAMIC_VISCOSITY values
prop = main_model_part.GetProperties()[1]
prop.SetValue(Kratos.DENSITY, 1)
prop.SetValue(Kratos.DYNAMIC_VISCOSITY, 0.002)

inletSubModelPart = main_model_part.CreateSubModelPart("Inlet")
outletSubModelPart = main_model_part.CreateSubModelPart("Outlet")

print(main_model_part.GetSubModelPart("LateralModelPart").NumberOfNodes()); input("PAUSE")
print(main_model_part.GetSubModelPart("LateralModelPart").NumberOfConditions()); input("PAUSE")

## Inlet and Outlet
for cond in main_model_part.GetSubModelPart("LateralModelPart").Conditions:
	print(cond.Id)
	nodes = cond.GetNodes()
	# UGLY! IMPROVE IT. We set the negative part as Inlet and the positive part as Outlet
	if (nodes[0].X <= 0 and nodes[1].X <= 0 and nodes[2].X <= 0):
		inletSubModelPart.AddNodes([nodes[0].Id, nodes[1].Id, nodes[2].Id])
		inletSubModelPart.AddConditions([cond.Id])
	else:
		outletSubModelPart.AddNodes([nodes[0].Id, nodes[1].Id, nodes[2].Id])
		outletSubModelPart.AddConditions([cond.Id])
print("*** inletSubModelPart and outletSubModelPart are filled ***")

print("\n*** ModelPart ***\n")

# we write mdpa file
mdpa_out_name = "cfd_data/test_{}/mdpa_file/11_Box_buildings_CFD".format(num_test)
building.WriteMdpaOutput(mdpa_out_name)
##################################################


""" DELETE IT """
# ##################################################
# # TODO: to add this part in geo_model.py

# # CFD PART
# main_model_part = building.GetGeoModelPart()

# current_model = Kratos.Model()

# if current_model.HasModelPart("NewModelPart"):
# 	# clear existing model part
# 	new_model_part = current_model.GetModelPart("NewModelPart")
# 	new_model_part.Elements.clear()
# 	new_model_part.Conditions.clear()
# 	new_model_part.Nodes.clear()
# else:
# 	new_model_part = current_model.CreateModelPart("NewModelPart")

# new_model_part.AddProperties(Kratos.Properties(1))

# # we set the DENSITY and DYNAMIC_VISCOSITY values
# prop = new_model_part.GetProperties()[1]
# prop.SetValue(Kratos.DENSITY, 1)
# prop.SetValue(Kratos.DYNAMIC_VISCOSITY, 0.002)

# fluidSubModelPart = new_model_part.CreateSubModelPart("Parts_Fluid")
# # wallSubModelPart = new_model_part.CreateSubModelPart("Wall")
# inletSubModelPart = new_model_part.CreateSubModelPart("Inlet")
# outletSubModelPart = new_model_part.CreateSubModelPart("Outlet")
# slipSubModelPart = new_model_part.CreateSubModelPart("Slip")
# print("*** The SubModelParts created ***")

# ## Fluid
# ### Nodes and Elements
# for elem in main_model_part.GetSubModelPart("ElementSubModelPart").Elements:
# 	nodes = elem.GetNodes()
# 	for node in nodes:
# 		n = new_model_part.CreateNewNode(node.Id, node.X, node.Y, node.Z)
# 		fluidSubModelPart.AddNode(n, 0)
# 	e = new_model_part.CreateNewElement("Element3D4N", elem.Id, [nodes[0].Id, nodes[1].Id, nodes[2].Id, nodes[3].Id], new_model_part.GetProperties()[1])
# 	fluidSubModelPart.AddElement(e, 0)
# print("*** Elements are created and are insert also in fluidSubModelPart ***")

# ## Inlet and Outlet
# for cond in main_model_part.GetSubModelPart("LateralModelPart").Conditions:
# 	nodes = cond.GetNodes()
# 	c = new_model_part.CreateNewCondition("WallCondition3D3N", cond.Id, [nodes[0].Id, nodes[1].Id, nodes[2].Id], new_model_part.GetProperties()[0])
# 	# UGLY! IMPROVE IT. We set the negative part as Inlet and the positive part as Outlet
# 	if (nodes[0].X <= 0 and nodes[1].X <= 0 and nodes[2].X <= 0):
# 		for node in nodes:
# 			n = new_model_part.CreateNewNode(node.Id, node.X, node.Y, node.Z)
# 			inletSubModelPart.AddNode(n, 0)
# 		inletSubModelPart.AddCondition(c, 0)
# 	else:
# 		for node in nodes:
# 			n = new_model_part.CreateNewNode(node.Id, node.X, node.Y, node.Z)
# 			outletSubModelPart.AddNode(n, 0)
# 		outletSubModelPart.AddCondition(c, 0)
# print("*** inletSubModelPart and outletSubModelPart are filled ***")

# ## Slip
# for cond in main_model_part.GetSubModelPart("BottomModelPart").Conditions:
# 	nodes = cond.GetNodes()
# 	for node in nodes:
# 		n = new_model_part.CreateNewNode(node.Id, node.X, node.Y, node.Z)
# 		slipSubModelPart.AddNode(n, 0)
	
# 	c = new_model_part.CreateNewCondition("WallCondition3D3N", cond.Id, [nodes[0].Id, nodes[1].Id, nodes[2].Id], new_model_part.GetProperties()[0])
# 	slipSubModelPart.AddCondition(c, 0)

# for cond in main_model_part.GetSubModelPart("TopModelPart").Conditions:
# 	nodes = cond.GetNodes()
# 	for node in nodes:
# 		n = new_model_part.CreateNewNode(node.Id, node.X, node.Y, node.Z)
# 		slipSubModelPart.AddNode(n, 0)
	
# 	c = new_model_part.CreateNewCondition("WallCondition3D3N", cond.Id, [nodes[0].Id, nodes[1].Id, nodes[2].Id], new_model_part.GetProperties()[0])
# 	slipSubModelPart.AddCondition(c, 0)
# print("*** slipSubModelPart is filled ***")

# # we write mdpa file
# mdpa_out_name = "cfd_data/test_{}/mdpa_file/11_Box_buildings_CFD".format(num_test)
# building.WriteMdpaOutput(mdpa_out_name)
# ##################################################


print("*** Time: ", time.time() - start_time)
print("\nTEST ", num_test, "END\n\n")
