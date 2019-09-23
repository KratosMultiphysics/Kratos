"""
# IT'S ONLY A TEST. SCRIPT NOT COMPLETE

	# edit:		06 September 2019
	# edit:		10 September 2019
	# edit:		12 September 2019 -> added import mdpa before step 09_...
	# edit:		20 September 2019 -> added GeoModel part
"""


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


def ComputeLimit(skin_model_part):
	
	Xmin = 1e7; Xmax = -1e7
	Ymin = 1e7; Ymax = -1e7
	Zmin = 1e7; Zmax = -1e7
	for node in skin_model_part.Nodes:
		if (node.X < Xmin): Xmin = node.X
		elif (node.X > Xmax): Xmax = node.X
		if (node.Y < Ymin): Ymin = node.Y
		elif (node.Y > Ymax): Ymax = node.Y
		if (node.Z < Zmin): Zmin = node.Z
		elif (node.Z > Zmax): Zmax = node.Z
	
	return (Xmin, Xmax, Ymin, Ymax, Zmin, Zmax)


start_time = time.time()
num_test = "55_september_simple_tetra_GeoModel"
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


importer = GeoImporter()
mesher = GeoMesher()
current_model = KratosMultiphysics.Model()

# skin model part (building)
skin_model_part = current_model.CreateModelPart("Skin")
skin_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISTANCE)
skin_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISTANCE_GRADIENT)
# nodes (building 1)
skin_model_part.CreateNewNode(1, 2.0, 2.0, -1e-7)
skin_model_part.CreateNewNode(2, 4.0, 2.0, -1e-7)
skin_model_part.CreateNewNode(3, 2.0, 4.0, -1e-7)
skin_model_part.CreateNewNode(4, 4.0, 4.0, -1e-7)
skin_model_part.CreateNewNode(5, 2.0, 2.0, 2.0)
skin_model_part.CreateNewNode(6, 4.0, 2.0, 2.0)
skin_model_part.CreateNewNode(7, 2.0, 4.0, 2.0)
skin_model_part.CreateNewNode(8, 4.0, 4.0, 2.0)
# nodes (building 2)
skin_model_part.CreateNewNode( 9, 6.0, 6.0, -1e-7)
skin_model_part.CreateNewNode(10, 8.0, 6.0, -1e-7)
skin_model_part.CreateNewNode(11, 6.0, 8.0, -1e-7)
skin_model_part.CreateNewNode(12, 8.0, 8.0, -1e-7)
skin_model_part.CreateNewNode(13, 6.0, 6.0, 2.0)
skin_model_part.CreateNewNode(14, 8.0, 6.0, 2.0)
skin_model_part.CreateNewNode(15, 6.0, 8.0, 2.0)
skin_model_part.CreateNewNode(16, 8.0, 8.0, 2.0)
# elements (building 1)
skin_model_part.CreateNewElement("Element2D3N",  1, [1, 5, 2], skin_model_part.GetProperties()[1])
skin_model_part.CreateNewElement("Element2D3N",  2, [2, 5, 6], skin_model_part.GetProperties()[1])
skin_model_part.CreateNewElement("Element2D3N",  3, [2, 6, 4], skin_model_part.GetProperties()[1])
skin_model_part.CreateNewElement("Element2D3N",  4, [4, 6, 8], skin_model_part.GetProperties()[1])
skin_model_part.CreateNewElement("Element2D3N",  5, [4, 8, 3], skin_model_part.GetProperties()[1])
skin_model_part.CreateNewElement("Element2D3N",  6, [3, 8, 7], skin_model_part.GetProperties()[1])
skin_model_part.CreateNewElement("Element2D3N",  7, [3, 7, 1], skin_model_part.GetProperties()[1])
skin_model_part.CreateNewElement("Element2D3N",  8, [1, 7, 5], skin_model_part.GetProperties()[1])
skin_model_part.CreateNewElement("Element2D3N",  9, [5, 7, 8], skin_model_part.GetProperties()[1])
skin_model_part.CreateNewElement("Element2D3N", 10, [5, 8, 6], skin_model_part.GetProperties()[1])
# elements (building 2)
skin_model_part.CreateNewElement("Element2D3N", 11, [ 9, 13, 10], skin_model_part.GetProperties()[1])
skin_model_part.CreateNewElement("Element2D3N", 12, [10, 13, 14], skin_model_part.GetProperties()[1])
skin_model_part.CreateNewElement("Element2D3N", 13, [10, 14, 12], skin_model_part.GetProperties()[1])
skin_model_part.CreateNewElement("Element2D3N", 14, [12, 14, 16], skin_model_part.GetProperties()[1])
skin_model_part.CreateNewElement("Element2D3N", 15, [12, 16, 11], skin_model_part.GetProperties()[1])
skin_model_part.CreateNewElement("Element2D3N", 16, [11, 16, 15], skin_model_part.GetProperties()[1])
skin_model_part.CreateNewElement("Element2D3N", 17, [11, 15,  9], skin_model_part.GetProperties()[1])
skin_model_part.CreateNewElement("Element2D3N", 18, [ 9, 15, 13], skin_model_part.GetProperties()[1])
skin_model_part.CreateNewElement("Element2D3N", 19, [13, 15, 16], skin_model_part.GetProperties()[1])
skin_model_part.CreateNewElement("Element2D3N", 20, [13, 16, 14], skin_model_part.GetProperties()[1])

# check if it is necessary
for node in skin_model_part.Nodes:
	node.Z -= 0.001

# save GiD file
mesher.SetGeoModelPart(skin_model_part)
mesher.CreateGidControlOutput("cfd_data/test_{}/gid_file/Ascii/03_skin_model_part".format(num_test), "GiD_PostAscii")
mesher.CreateGidControlOutput("cfd_data/test_{}/gid_file/Binary/03_skin_model_part".format(num_test), "GiD_PostBinary")
# writing file mdpa
mdpa_out_name = "cfd_data/test_{}/mdpa_file/03_skin_model_part".format(num_test)
mesher.WriteMdpaOutput(mdpa_out_name)

# compute the limit of the domain
Xmin, Xmax, Ymin, Ymax, Zmin, Zmax = ComputeLimit(skin_model_part)

incr = 10	# domain offset
# new function (under test)
importer.CreateBoxTetra("Volume", Xmin-incr, Xmax+incr, Ymin-incr, Ymax+incr, 0.0, Zmax+incr)
model_part = importer.GetGeoModelPart()
# save GiD file
mesher.SetGeoModelPart(model_part)
mesher.CreateGidControlOutput("cfd_data/test_{}/gid_file/Ascii/01_box".format(num_test), "GiD_PostAscii")
mesher.CreateGidControlOutput("cfd_data/test_{}/gid_file/Binary/01_box".format(num_test), "GiD_PostBinary")
# writing file mdpa
mdpa_out_name = "cfd_data/test_{}/mdpa_file/01_box".format(num_test)
mesher.WriteMdpaOutput(mdpa_out_name)

# refine mesh with MMG
mesher.SetGeoModelPart(model_part)
mesher.RefineMesh_test(1.0, 10.0, 1.0)
mesher.CreateGidControlOutput("cfd_data/test_{}/gid_file/Ascii/02_box_refinement_1".format(num_test), "GiD_PostAscii")
mesher.CreateGidControlOutput("cfd_data/test_{}/gid_file/Binary/02_box_refinement_1".format(num_test), "GiD_PostBinary")
# writing file mdpa
mdpa_out_name = "cfd_data/test_{}/mdpa_file/02_box_refinement_1".format(num_test)
mesher.WriteMdpaOutput(mdpa_out_name)


# building steps
print("\n\nSTART PROCESS BUILDING")
stop_begin_loop = time.time()

building = GeoBuilding()
building.SetGeoModelPart(model_part)
building.ImportBuilding(skin_model_part)

# # 1st cut
# print("\n\n***** STEP 1 *****\n")
# # distance field from hull
# building.ComputeDistanceFieldFromHull(False, 1e-7)
# building.CreateGidControlOutput("cfd_data/test_{}/gid_file/Ascii/04_Box_buildings_distance_field_1".format(num_test), "GiD_PostAscii")
# building.CreateGidControlOutput("cfd_data/test_{}/gid_file/Binary/04_Box_buildings_distance_field_1".format(num_test), "GiD_PostBinary")
# # writing file mdpa
# mdpa_out_name = "cfd_data/test_{}/mdpa_file/04_Box_buildings_distance_field_1".format(num_test)
# building.WriteMdpaOutput(mdpa_out_name)
# # subtract buildings
# building.SubtractBuildingMOD(0.5, 10.0, 0.1, "Linear", "STANDARD", "false")			# interpolation = constant; disc_type = STANDARD; remove_regions=false
# building.CreateGidControlOutput("cfd_data/test_{}/gid_file/Ascii/05_Box_buildings_subtracted_1".format(num_test), "GiD_PostAscii")
# building.CreateGidControlOutput("cfd_data/test_{}/gid_file/Binary/05_Box_buildings_subtracted_1".format(num_test), "GiD_PostBinary")
# # writing file mdpa
# mdpa_out_name = "cfd_data/test_{}/mdpa_file/05_Box_buildings_subtracted_1".format(num_test)
# building.WriteMdpaOutput(mdpa_out_name)

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
# # building.SubtractBuildingMOD(0.2, 10.0, 0.1, "Constant", "STANDARD", "false")			# interpolation = constant; disc_type = STANDARD; remove_regions=false
# building.SubtractBuildingMOD(0.1, 10.0, 0.1, "Constant", "STANDARD", "false")			# interpolation = constant; disc_type = STANDARD; remove_regions=false
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
# subtract buildings
# building.SubtractBuildingMOD(0.2, 10.0, 0.1, "exponential", "ISOSURFACE", "true")			# interpolation = exponential; disc_type = STANDARD; remove_regions=true
building.SubtractBuildingMOD(0.1, 10.0, 0.1, "exponential", "ISOSURFACE", "true")			# interpolation = exponential; disc_type = STANDARD; remove_regions=true

# we write the GiD file
building.CreateGidControlOutput("cfd_data/test_{}/gid_file/Ascii/09_Box_buildings_subtracted_3_1".format(num_test), "GiD_PostAscii")
building.CreateGidControlOutput("cfd_data/test_{}/gid_file/Binary/09_Box_buildings_subtracted_3_1".format(num_test), "GiD_PostBinary")

stop_01 = time.time()
# we delete the condition if at least one node it is not in main model part
main_model_part = building.GetGeoModelPart()
KratosGeo.CleaningUtilities(main_model_part).CleanConditions()
stop_02 = time.time()
print("CleanConditions done in {} s".format(stop_02-stop_01))

# we set the geo model part to save the GiD file and mdpa file
building.SetGeoModelPart(main_model_part)

# we write the GiD file
building.CreateGidControlOutput("cfd_data/test_{}/gid_file/Ascii/10_Box_buildings_after_CleanConditions".format(num_test), "GiD_PostAscii")
building.CreateGidControlOutput("cfd_data/test_{}/gid_file/Binary/10_Box_buildings_after_CleanConditions".format(num_test), "GiD_PostBinary")
# writing file mdpa
mdpa_out_name = "cfd_data/test_{}/mdpa_file/10_Box_buildings_after_CleanConditions".format(num_test)
building.WriteMdpaOutput(mdpa_out_name)
print("*************************************************************************")

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


# we copy MainKratos.py file
from shutil import copyfile
src_path = "data/MainKratos.py"
dest_path = "cfd_data/test_{}/analysis_file/MainKratos.py".format(num_test)
copyfile(src_path, dest_path)

# we copy json file parameters
src_path = "data/parameters/ProjectParameters.json"
dest_path = "cfd_data/test_{}/analysis_file/ProjectParameters.json".format(num_test)
copyfile(src_path, dest_path)


#
lateral_model_part = building.GetGeoModelPart().GetSubModelPart("LateralModelPart")
inlet_model_part = building.GetGeoModelPart().CreateSubModelPart("Inlet")
outlet_model_part = building.GetGeoModelPart().CreateSubModelPart("Outlet")
slip_model_part = building.GetGeoModelPart().CreateSubModelPart("Slip")
for cond in lateral_model_part.Conditions:
	nodes = cond.GetNodes()
	if (nodes[0].X == -8.0) and (nodes[1].X == -8.0) and (nodes[2].X == -8.0):
		inlet_model_part.AddCondition(cond, 0)
		for node in nodes:
			inlet_model_part.AddNode(node, 0)
	elif (nodes[0].X > 17.9) and (nodes[1].X > 17.9) and (nodes[2].X > 17.9):
		outlet_model_part.AddCondition(cond, 0)
		for node in nodes:
			outlet_model_part.AddNode(node, 0)
	else:
		slip_model_part.AddCondition(cond, 0)
		for node in nodes:
			slip_model_part.AddNode(node, 0)

building.GetGeoModelPart().RemoveSubModelPart("LateralModelPart")

print("\n*** ModelPart ***\n", building.GetGeoModelPart())
#

print("main_model_part\n", main_model_part)

model = GeoModel()
model.SetGeoModelPart(main_model_part)
model.GenerateCfdModelPart()

# print("\nGeoCfdModelPart BEFORE")
# print(model.GetGeoCfdModelPart())
print("*****************")

stop_05 = time.time()
model.FillPartsFluid("Parts_Fluid")
print("Parts_Fluid DONE")

model.FillNoslip("SKIN_ISOSURFACE")
print("SKIN_ISOSURFACE DONE")

model.FillSlip("TopModelPart")
print("TopModelPart DONE")

model.FillSlip("BottomModelPart")
print("BottomModelPart DONE")

model.FillSlip("Slip")
print("Slip DONE")

model.FillInlet("Inlet")
print("Inlet DONE")

model.FillOutlet("Outlet")
print("Outlet DONE")

# ##################################### call C++ function
# stop_03 = time.time()
# # FillCfdModelpartUtilities
# CfdModelPart = model.GetGeoCfdModelPart()
# KratosGeo.FillCfdModelpartUtilities(CfdModelPart).FillPartsFluid(main_model_part, "Parts_Fluid")
# KratosGeo.FillCfdModelpartUtilities(CfdModelPart).FillNoslip(main_model_part, "SKIN_ISOSURFACE")
# KratosGeo.FillCfdModelpartUtilities(CfdModelPart).FillSlip(main_model_part, "TopModelPart")
# KratosGeo.FillCfdModelpartUtilities(CfdModelPart).FillSlip(main_model_part, "BottomModelPart")
# KratosGeo.FillCfdModelpartUtilities(CfdModelPart).FillSlip(main_model_part, "Slip")
# KratosGeo.FillCfdModelpartUtilities(CfdModelPart).FillInlet(main_model_part, "Inlet")
# KratosGeo.FillCfdModelpartUtilities(CfdModelPart).FillOutlet(main_model_part, "Outlet")
# stop_04 = time.time()
# print("FillPartsFluid done in {} s".format(stop_04-stop_03))

# print("\n*****************************\n", model.GetGeoCfdModelPart(), "\n*****************************\n")
# ##################################### end call C++ function


stop_06 = time.time()
print("\n* Filled time: ", stop_06-stop_05)

print("\nGeoCfdModelPart AFTER")
print(model.GetGeoCfdModelPart())

# writing file mdpa
model.SetGeoModelPart(model.GetGeoCfdModelPart())
mdpa_out_name = "cfd_data/test_{}/analysis_file/12_CFD_model".format(num_test)
model.WriteMdpaOutput(mdpa_out_name)

print("*** Time: ", time.time() - start_time)
print("\nTEST ", num_test, "END\n\n")
