"""
	# create:	09 September 2019 -> we import the mdpa file and we fill a new model part to perform a CFD analysis
	# edit:		10 September 2019
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
num_test = "44_september_import_and_fill_CFD_modelpart"
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


# import domain from mdpa file
importer._InitializeModelPart("test_model")
terrain_model_part = importer.ModelPart

model_part_in = "data/mdpa_file/domain_13_sept_2019/09_Box_buildings_subtracted_3_after_CleanConditions"
KratosMultiphysics.ModelPartIO(model_part_in).ReadModelPart(terrain_model_part)
print("\n*** MDPA imported! ***\n")

building.SetGeoModelPart(terrain_model_part)
building.CreateGidControlOutput("cfd_data/test_{}/gid_file/Ascii/09_Box_buildings_subtracted_3_after_CleanConditions".format(num_test), "GiD_PostAscii")
building.CreateGidControlOutput("cfd_data/test_{}/gid_file/Binary/09_Box_buildings_subtracted_3_after_CleanConditions".format(num_test), "GiD_PostBinary")

main_model_part = building.GetGeoModelPart()
inletSubModelPart = main_model_part.CreateSubModelPart("Inlet")
outletSubModelPart = main_model_part.CreateSubModelPart("Outlet")
slipSubModelPart = main_model_part.CreateSubModelPart("Slip")

inlet_nodes = [];	inlet_cond = []
outlet_nodes = [];	outlet_cond =[]
slip_nodes = [];	slip_cond = []
for cond in main_model_part.GetSubModelPart("LateralModelPart").Conditions:
	nodes = cond.GetNodes()
	if (nodes[0].X == -8.0 and nodes[1].X == -8.0 and nodes[2].X == -8.0):
		inlet_nodes.extend([nodes[0].Id, nodes[1].Id, nodes[2].Id])
		inlet_cond.append(cond.Id)

	elif (nodes[0].X > 17.9 and nodes[1].X > 17.9 and nodes[2].X > 17.9):
		outlet_nodes.extend([nodes[0].Id, nodes[1].Id, nodes[2].Id])
		outlet_cond.append(cond.Id)
	
	else:
		slip_nodes.extend([nodes[0].Id, nodes[1].Id, nodes[2].Id])
		slip_cond.append(cond.Id)

inletSubModelPart.AddNodes(inlet_nodes)
inletSubModelPart.AddConditions(inlet_cond)

outletSubModelPart.AddNodes(outlet_nodes)
outletSubModelPart.AddConditions(outlet_cond)

slipSubModelPart.AddNodes(slip_nodes)
slipSubModelPart.AddConditions(slip_cond)


print(main_model_part)
building.SetGeoModelPart(main_model_part)
# writing file mdpa
mdpa_out_name = "cfd_data/test_{}/mdpa_file/10_Box_buildings".format(num_test)
building.WriteMdpaOutput(mdpa_out_name)



# ##################################################
# # CFD PART
# # TODO: to add this part in geo_model.py
# main_model_part = building.GetGeoModelPart()


### TOO SLOW ###

# current_model = KratosMultiphysics.Model()

# if current_model.HasModelPart("NewModelPart"):
# 	# clear existing model part
# 	new_model_part = current_model.GetModelPart("NewModelPart")
# 	new_model_part.Elements.clear()
# 	new_model_part.Conditions.clear()
# 	new_model_part.Nodes.clear()
# else:
# 	new_model_part = current_model.CreateModelPart("NewModelPart")

# new_model_part.AddProperties(KratosMultiphysics.Properties(1))

# # we set the DENSITY and DYNAMIC_VISCOSITY values
# prop = new_model_part.GetProperties()[1]
# prop.SetValue(KratosMultiphysics.DENSITY, 1)
# prop.SetValue(KratosMultiphysics.DYNAMIC_VISCOSITY, 0.002)

# fluidSubModelPart = new_model_part.CreateSubModelPart("Parts_Fluid")
# # wallSubModelPart = new_model_part.CreateSubModelPart("Wall")
# inletSubModelPart = new_model_part.CreateSubModelPart("Inlet")
# outletSubModelPart = new_model_part.CreateSubModelPart("Outlet")
# slipSubModelPart = new_model_part.CreateSubModelPart("Slip")
# print("*** The SubModelParts created ***")

# ## Fluid
# ### Nodes and Elements
# time_01 = time.time()
# node_list = []
# elem_list = []
# for elem in main_model_part.GetSubModelPart("ElementSubModelPart").Elements:
# 	nodes = elem.GetNodes()
# 	for node in nodes:
# 		n = new_model_part.CreateNewNode(node.Id, node.X, node.Y, node.Z)
# 		# fluidSubModelPart.AddNode(n, 0)
# 		node_list.append(n.Id)
# 	e = new_model_part.CreateNewElement("Element3D4N", elem.Id, [nodes[0].Id, nodes[1].Id, nodes[2].Id, nodes[3].Id], new_model_part.GetProperties()[1])
# 	# fluidSubModelPart.AddElement(e, 0)
# 	elem_list.append(e.Id)
# fluidSubModelPart.AddNodes(node_list)
# fluidSubModelPart.AddElements(elem_list)
# time_02 = time.time()
# print("*** Elements are created and are insert also in fluidSubModelPart in: {} ***".format(time_02 - time_01))

# ## Inlet and Outlet
# time_03 = time.time()
# node_inlet_list = []
# cond_inlet_list = []
# node_outlet_list = []
# cond_outlet_list = []
# for cond in main_model_part.GetSubModelPart("LateralModelPart").Conditions:
# 	nodes = cond.GetNodes()
# 	c = new_model_part.CreateNewCondition("WallCondition3D3N", cond.Id, [nodes[0].Id, nodes[1].Id, nodes[2].Id], new_model_part.GetProperties()[0])
# 	# UGLY! IMPROVE IT. We set the negative part as Inlet and the positive part as Outlet
# 	if (nodes[0].X <= 0 and nodes[1].X <= 0 and nodes[2].X <= 0):
# 		for node in nodes:
# 			n = new_model_part.CreateNewNode(node.Id, node.X, node.Y, node.Z)
# 			# inletSubModelPart.AddNode(n, 0)
# 			node_inlet_list.append(n.Id)
# 		# inletSubModelPart.AddCondition(c, 0)
# 		cond_inlet_list.append(c.Id)
# 	else:
# 		for node in nodes:
# 			n = new_model_part.CreateNewNode(node.Id, node.X, node.Y, node.Z)
# 			# outletSubModelPart.AddNode(n, 0)
# 			node_outlet_list.append(n.Id)
# 		# outletSubModelPart.AddCondition(c, 0)
# 		cond_outlet_list.append(c.Id)
# inletSubModelPart.AddNodes(node_inlet_list)
# inletSubModelPart.AddConditions(cond_inlet_list)
# outletSubModelPart.AddNodes(node_outlet_list)
# outletSubModelPart.AddConditions(cond_outlet_list)
# time_04 = time.time()
# print("*** inletSubModelPart and outletSubModelPart are filled in: {} ***".format(time_04 - time_03))

# ## Slip
# time_05 = time.time()
# node_list = []
# cond_list = []
# for cond in main_model_part.GetSubModelPart("BottomModelPart").Conditions:
# 	nodes = cond.GetNodes()
# 	for node in nodes:
# 		n = new_model_part.CreateNewNode(node.Id, node.X, node.Y, node.Z)
# 		# slipSubModelPart.AddNode(n, 0)
# 		node_list.append(n.Id)
	
# 	c = new_model_part.CreateNewCondition("WallCondition3D3N", cond.Id, [nodes[0].Id, nodes[1].Id, nodes[2].Id], new_model_part.GetProperties()[0])
# 	# slipSubModelPart.AddCondition(c, 0)
# 	cond_list.append(c.Id)

# for cond in main_model_part.GetSubModelPart("TopModelPart").Conditions:
# 	nodes = cond.GetNodes()
# 	for node in nodes:
# 		n = new_model_part.CreateNewNode(node.Id, node.X, node.Y, node.Z)
# 		# slipSubModelPart.AddNode(n, 0)
# 		node_list.append(n.Id)
	
# 	c = new_model_part.CreateNewCondition("WallCondition3D3N", cond.Id, [nodes[0].Id, nodes[1].Id, nodes[2].Id], new_model_part.GetProperties()[0])
# 	# slipSubModelPart.AddCondition(c, 0)
# 	cond_list.append(c.Id)
# slipSubModelPart.AddNodes(node_list)
# slipSubModelPart.AddConditions(cond_list)
# time_06 = time.time()
# print("*** slipSubModelPart is filled in: {} ***".format(time_06 - time_05))

# # we write mdpa file
# mdpa_out_name = "cfd_data/test_{}/mdpa_file/11_Box_buildings_CFD".format(num_test)
# building.WriteMdpaOutput(mdpa_out_name)
# ##################################################


print("*** Time: ", time.time() - start_time)
print("\nTEST ", num_test, "END\n\n")

