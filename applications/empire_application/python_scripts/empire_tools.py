from __future__ import print_function, absolute_import, division # makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# import libraries
from KratosMultiphysics import *
import ctypes as ctp
import os
import numpy as empire_numpy


class EmpireTools:
    
    def __init__(self,model_part):
      self.model_part = model_part
      self.number_of_nodes = self.model_part.NumberOfNodes()
      
      
    def ExtractInterface(self,numNodes, numElems, nodeCoors, nodeIDs, numNodesPerElem, elemTable):
      numNodes.append(self.number_of_nodes)
      numElems.append(self.model_part.NumberOfConditions())
      
      for node in self.model_part.Nodes:
        nodeCoors.append(node.X)
        nodeCoors.append(node.Y)
        nodeCoors.append(node.Z)
        nodeIDs.append(node.Id)
        
      for cond in self.model_part.Conditions:
        numNodesPerElem.append(len(cond.GetNodes()))
        for node in cond.GetNodes():
            elemTable.append(node.Id)
          
    def GetDataField(self, data_field_name, data_field):
        for node in self.model_part.Nodes:
            data_value = node.GetSolutionStepValue(data_field_name,0)
            for i in range(len(data_value)):
                data_field.append(data_value[i])

    def SetDataField(self, data_field_name, data_field):
        # TODO somehow get the size of the variable
        NodesArray

        # check if size of data field is correct
        if len(data_field) != self.number_of_nodes * size_of_variable:
            raise("ERROR: received data field has wrong size!")

        value = Vector(size_of_variable)
        i = 0
        # assign values to nodes of interface for current time step
        for node in self.model_part.Nodes:
            for j in range(size_of_variable):
                value[j] = data_field[size_of_variable * i + j]
            
            node.SetSolutionStepValue(kratos_variable, 0, value)
            
            i = i + 1