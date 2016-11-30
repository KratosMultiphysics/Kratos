from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

from KratosMultiphysics import *
from KratosMultiphysics.MappingApplication import *

import os
import process_factory
import KratosMultiphysics.KratosUnittest as KratosUnittest

class KratosExecuteNearestNeighborMapperTest(KratosUnittest.TestCase):

    def __init__(self):
        # Mdpa Input files
        input_file_structure = "Mapper_Test_1/FSI_Example4Mapper_1_Structural"
        input_file_fluid     = "Mapper_Test_1/FSI_Example4Mapper_1_Fluid"

        variable_list = [PRESSURE, VELOCITY]
        self.model_part_origin  = self.ReadModelPartSerial("ModelPartNameOrigin", input_file_fluid, variable_list)
        self.model_part_destination = self.ReadModelPartSerial("ModelPartNameDestination", input_file_structure, variable_list)

        self.interface_sub_model_part_origin = self.model_part_origin.GetSubModelPart("FluidNoSlipInterface3D_interface_orig_fluid")
        self.interface_sub_model_part_destination = self.model_part_destination.GetSubModelPart("StructureInterface3D_interface_dest_struct")

        # Initialize Mapper
        self.nearestNeighborMapper = NearestNeighborMapper(self.interface_sub_model_part_origin, self.interface_sub_model_part_destination)


    def TestMapConstantScalarValues(self):
        map_value = 5.123
        variable_origin = PRESSURE
        variable_destination = PRESSURE

        self.SetValuesOnNodes(self.model_part_origin, variable_origin, map_value)

        # Overwriting Values
        self.nearestNeighborMapper.Map(variable_origin, variable_destination)
        self.CheckValues(self.interface_sub_model_part_destination, variable_destination, map_value)

        # Adding Values
        self.nearestNeighborMapper.Map(variable_origin, variable_destination, True)
        self.CheckValues(self.interface_sub_model_part_destination, variable_destination, map_value*2)

    def TestInverseMapConstantScalarValues(self):
        map_value = -8.6647
        variable_origin = PRESSURE
        variable_destination = PRESSURE

        self.SetValuesOnNodes(self.model_part_destination, variable_destination, map_value)

        # Overwriting Values
        self.nearestNeighborMapper.InverseMap(variable_origin, variable_destination)
        self.CheckValues(self.interface_sub_model_part_origin, variable_origin, map_value)

        # Adding Values
        self.nearestNeighborMapper.InverseMap(variable_origin, variable_destination, True)
        self.CheckValues(self.interface_sub_model_part_origin, variable_origin, map_value*2)

    def TestMapConstantVectorValues(self):
        map_value = [15.99, -2.88, 3.123]
        variable_origin = VELOCITY
        variable_destination = VELOCITY

        self.SetValuesOnNodes(self.model_part_origin, VELOCITY_X, map_value[0])
        self.SetValuesOnNodes(self.model_part_origin, VELOCITY_Y, map_value[1])
        self.SetValuesOnNodes(self.model_part_origin, VELOCITY_Z, map_value[2])

        # Overwriting Values
        self.nearestNeighborMapper.Map(variable_origin, variable_destination)

        self.CheckValues(self.interface_sub_model_part_destination, VELOCITY_X, map_value[0])
        self.CheckValues(self.interface_sub_model_part_destination, VELOCITY_Y, map_value[1])
        self.CheckValues(self.interface_sub_model_part_destination, VELOCITY_Z, map_value[2])

        # Adding Values
        self.nearestNeighborMapper.Map(variable_origin, variable_destination, True)

        self.CheckValues(self.interface_sub_model_part_destination, VELOCITY_X, map_value[0]*2)
        self.CheckValues(self.interface_sub_model_part_destination, VELOCITY_Y, map_value[1]*2)
        self.CheckValues(self.interface_sub_model_part_destination, VELOCITY_Z, map_value[2]*2)

    def TestInverseMapConstantVectorValues(self):
        map_value = [1.4785, -0.88, -33.123]
        variable_origin = VELOCITY
        variable_destination = VELOCITY

        self.SetValuesOnNodes(self.model_part_destination, VELOCITY_X, map_value[0])
        self.SetValuesOnNodes(self.model_part_destination, VELOCITY_Y, map_value[1])
        self.SetValuesOnNodes(self.model_part_destination, VELOCITY_Z, map_value[2])

        # Overwriting Values
        self.nearestNeighborMapper.InverseMap(variable_origin, variable_destination)

        self.CheckValues(self.interface_sub_model_part_origin, VELOCITY_X, map_value[0])
        self.CheckValues(self.interface_sub_model_part_origin, VELOCITY_Y, map_value[1])
        self.CheckValues(self.interface_sub_model_part_origin, VELOCITY_Z, map_value[2])

        # Adding Values
        self.nearestNeighborMapper.InverseMap(variable_origin, variable_destination, True)

        self.CheckValues(self.interface_sub_model_part_origin, VELOCITY_X, map_value[0]*2)
        self.CheckValues(self.interface_sub_model_part_origin, VELOCITY_Y, map_value[1]*2)
        self.CheckValues(self.interface_sub_model_part_origin, VELOCITY_Z, map_value[2]*2)


    def TestMapNonConstantScalarValues(self):
        variable_origin = PRESSURE
        variable_destination = PRESSURE

        self.SetValuesOnNodesId(self.interface_sub_model_part_origin, variable_origin)

        self.nearestNeighborMapper.Map(variable_origin, variable_destination)
        nodal_values = [12,14,13,22,25,23,42,59,60,57,76,77,75,89,93,112,113,115,114,121,124,123,148,147,154]
        self.CheckValuesId(self.interface_sub_model_part_destination, variable_destination, nodal_values)

    def TestInverseMapNonConstantScalarValues(self):
        variable_origin = PRESSURE
        variable_destination = PRESSURE

        self.SetValuesOnNodesId(self.interface_sub_model_part_destination, variable_destination)

        self.nearestNeighborMapper.InverseMap(variable_origin, variable_destination)
        nodal_values = [11,17,16,26,33,35,32,50,53,35,33,32,50,54,50,53,74,66,67,86,117,86,96,96,110,111,118,117,120,96,132,131,132,131,146,145,153]
        self.CheckValuesId(self.interface_sub_model_part_origin, variable_origin, nodal_values)


    def ReadModelPartSerial(self, model_part_name, model_part_input_file, variable_list):
        model_part = ModelPart(model_part_name)
        for variable in variable_list:
            model_part.AddNodalSolutionStepVariable(variable)

        model_part_io = ModelPartIO(model_part_input_file)
        model_part_io.ReadModelPart(model_part)

        model_part.ProcessInfo.SetValue(DOMAIN_SIZE, 3)
        model_part.SetBufferSize(1)

        return model_part

    def SetValuesOnNodes(self, model_part, variable, value):
        for node in model_part.Nodes:
            node.SetSolutionStepValue(variable, value)

    def SetValuesOnNodesId(self, model_part, variable):
        for node in model_part.Nodes:
            node.SetSolutionStepValue(variable, node.Id)

    def CheckValues(self, model_part, variable, map_value):
        for node in model_part.Nodes:
            value = node.GetSolutionStepValue(variable)
            self.assertAlmostEqual(map_value,value)

    def CheckValuesId(self, model_part, variable, nodal_values):
        if (len(nodal_values) != len(model_part.Nodes)):
            raise Exception("Array Sizes are not matching")

        i = 0
        for node in model_part.Nodes:
            value = node.GetSolutionStepValue(variable)
            self.assertAlmostEqual(nodal_values[i],value)
            i += 1
