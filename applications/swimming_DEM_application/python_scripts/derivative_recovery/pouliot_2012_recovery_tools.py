from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
import math
import os
from KratosMultiphysics import *
from KratosMultiphysics.SwimmingDEMApplication import *

class PouliotRecoveryTools:
    def __init__(self):
        pass

    def MakeRecoveryModelPart(self, model_part):
        edges_model_part = ModelPart("Edges")
        set_of_all_edges = set()
        for elem in model_part.Elements:
            for i, first_node in enumerate(elem.Nodes[:-1]):
                for j, second_node in enumerate(elem.Nodes[i:]):
                    edge_ids = (first_node.Id, second_node.Id)
                    set_of_all_edges.add(edge_ids)
        for i, edge in enumerate(set_of_all_edges):
            edges_model_part.CreateNewElement("Element3D1N", i, edge, edges_model_part.GetProperties()[0])
