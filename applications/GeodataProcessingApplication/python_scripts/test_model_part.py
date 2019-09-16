import KratosMultiphysics as Kratos
import KratosMultiphysics.GeodataProcessingApplication as KratosGeo
import KratosMultiphysics.FluidDynamicsApplication

from geo_importer import GeoImporter
from geo_mesher import GeoMesher
from geo_preprocessor import GeoPreprocessor
from geo_building import GeoBuilding

import time
import os


current_model = KratosMultiphysics.Model()
ModelPart = current_model.CreateModelPart("name_model_part")

SMP_1 = ModelPart.CreateSubModelPart("SMP_1")

node1 = ModelPart.CreateNewNode(1, 0.0, 0.0, 0.0)
node2 = ModelPart.CreateNewNode(2, 1.0, 0.0, 0.0)
node3 = ModelPart.CreateNewNode(3, 1.0, 1.0, 0.0)
node4 = ModelPart.CreateNewNode(4, 0.0, 0.0, 1.0)

elem = ModelPart.CreateNewElement("Element3D4N",  1, [1, 2, 3, 4], ModelPart.GetProperties()[1])

# print(dir(elem))
print(elem.GetGeometry())

# SMP_1.AddNode(node, 0)

# for sub_model in ModelPart.SubModelParts:
# 	for node in sub_model.Nodes:
# 		node.Set(KratosMultiphysics.TO_ERASE, True)

# # we erase unused nodes
# ModelPart.RemoveNodesFromAllLevels(KratosMultiphysics.TO_ERASE)

# print(ModelPart)