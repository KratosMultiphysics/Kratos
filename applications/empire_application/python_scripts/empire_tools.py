from __future__ import print_function, absolute_import, division # makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# import libraries
from KratosMultiphysics import *
import ctypes as ctp
import os
import numpy as empire_numpy


class EmpireTools:
    
    def __init__(self,model_part):
      self.model_part = model_part      
      
    def ExtractInterface(self,num_nodes, num_elements, node_coords, node_IDs, num_nodes_per_element, element_table):
      num_nodes.append(self.model_part.NumberOfNodes())
      num_elements.append(self.model_part.NumberOfConditions())
      
      for node in self.model_part.Nodes:
        node_coords.append(node.X)
        node_coords.append(node.Y)
        node_coords.append(node.Z)
        node_IDs.append(node.Id)
        
      for cond in self.model_part.Conditions:
        num_nodes_per_element.append(len(cond.GetNodes()))
        for node in cond.GetNodes():
            element_table.append(node.Id)

    def ConstructMesh(self, num_nodes, num_elements, node_coords, node_IDs, num_nodes_per_element, element_table, construct_conditions):
        # This function requires an empty ModelPart
        # It constructs Nodes and Conditions from what was received from Empire

        # Some checks to validate input:
        if self.model_part.NumberOfNodes() != 0:
            raise Exception("ModelPart is not empty, it has some Nodes!")

        # Create Nodes
        for i in range(num_nodes):
            self.model_part.CreateNewNode(node_IDs[i], node_coords[3*i+0], node_coords[3*i+1], node_coords[3*i+2]) # Id, X, Y, Z

        if construct_conditions:
            # Create Property for Condition
            self.model_part.AddProperties(Properties(1))
            prop = self.model_part.GetProperties()[1]

            element_table_counter = 0
            # Create Conditions
            for i in range(num_elements):
                num_nodes_condition = num_nodes_per_element[i]
                if num_nodes_condition == 2:
                    name_condition = "LineCondition3D2N"
                elif num_nodes_condition == 3:
                    name_condition = "SurfaceCondition3D3N"
                elif num_nodes_condition == 4:
                    name_condition = "SurfaceCondition3D4N"
                else:
                    raise Exception("Wrong number of nodes for creating the condition")

                condition_nodes = []
                for j in range(num_nodes_condition):
                    condition_nodes.append(int(element_table[element_table_counter]))
                    element_table_counter += 1
                  
                self.model_part.CreateNewCondition(name_condition, i+1, condition_nodes, prop)

    def GetDataField(self, data_field_name, data_field):
        for node in self.model_part.Nodes:
            data_value = node.GetSolutionStepValue(data_field_name)
            for i in range(len(data_value)):
                data_field.append(data_value[i])

    def SetDataField(self, data_field_name, data_field, size_of_variable):
        # check if size of data field is correct
        if len(data_field) != self.model_part.NumberOfNodes() * size_of_variable:
            raise("ERROR: received data field has wrong size!")

        value = Vector(size_of_variable)
        i = 0
        # assign values to nodes of interface for current time step
        for node in self.model_part.Nodes:
            for j in range(size_of_variable):
                value[j] = data_field[size_of_variable * i + j]
            
            node.SetSolutionStepValue(data_field_name, 0, value)
            
            i = i + 1