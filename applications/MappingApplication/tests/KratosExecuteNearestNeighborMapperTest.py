from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

from KratosMultiphysics import *
from KratosMultiphysics.MappingApplication import *
import Mapper

import os
import process_factory
import KratosMultiphysics.KratosUnittest as KratosUnittest

class KratosExecuteNearestNeighborMapperTest(KratosUnittest.TestCase):

    def __init__(self, gid_output):
        self.GiD_output = gid_output

        # Mdpa Input files
        input_file_structure = "KratosExecuteNearestNeighborMapperTest_mdpa/FSI_Example4Mapper_1_Structural"
        input_file_fluid     = "KratosExecuteNearestNeighborMapperTest_mdpa/FSI_Example4Mapper_1_Fluid"

        parameter_file = open("KratosExecuteNearestNeighborMapperTest_mdpa/KratosExecuteNearestNeighborMapperTest.json",'r')
        ProjectParameters = Parameters(parameter_file.read())

        variable_list = [PRESSURE, VELOCITY]
        self.model_part_origin  = self.ReadModelPartSerial("ModelPartNameOrigin", input_file_fluid, variable_list)
        self.model_part_destination = self.ReadModelPartSerial("ModelPartNameDestination", input_file_structure, variable_list)

        # needed for the tests only, usually one does not need to get the submodel-parts for the mapper
        self.interface_sub_model_part_origin = self.model_part_origin.GetSubModelPart("FluidNoSlipInterface3D_interface_orig_fluid")
        self.interface_sub_model_part_destination = self.model_part_destination.GetSubModelPart("StructureInterface3D_interface_dest_struct")

        # Initialize Mapper
        self.nearestNeighborMapper = Mapper.NonMatchingGridMapper(self.model_part_origin, self.model_part_destination, ProjectParameters)

        if (self.GiD_output):
            self.InitializeGiD()


    def TestMapConstantScalarValues(self, output_time):
        map_value = 5.123
        variable_origin = PRESSURE
        variable_destination = PRESSURE

        self.SetValuesOnNodes(self.model_part_origin, variable_origin, map_value)

        if (self.GiD_output):
            self.WriteNodalResultsCustom(self.gid_io_origin, self.model_part_origin, variable_origin, output_time)

        # Overwriting Values
        self.nearestNeighborMapper.Map(variable_origin, variable_destination)

        if (self.GiD_output):
            self.WriteNodalResultsCustom(self.gid_io_destination, self.model_part_destination, variable_destination, output_time)

        self.CheckValues(self.interface_sub_model_part_destination, variable_destination, map_value)

        # Adding Values
        self.nearestNeighborMapper.Map(variable_origin, variable_destination, True)

        if (self.GiD_output):
            self.WriteNodalResultsCustom(self.gid_io_destination, self.model_part_destination, variable_destination, output_time + 0.1)

        self.CheckValues(self.interface_sub_model_part_destination, variable_destination, map_value*2)

    def TestInverseMapConstantScalarValues(self, output_time):
        map_value = -8.6647
        variable_origin = PRESSURE
        variable_destination = PRESSURE

        self.SetValuesOnNodes(self.model_part_destination, variable_destination, map_value)

        if (self.GiD_output):
            self.WriteNodalResultsCustom(self.gid_io_destination, self.model_part_destination, variable_destination, output_time)

        # Overwriting Values
        self.nearestNeighborMapper.InverseMap(variable_origin, variable_destination)

        if (self.GiD_output):
            self.WriteNodalResultsCustom(self.gid_io_origin, self.model_part_origin, variable_origin, output_time)

        self.CheckValues(self.interface_sub_model_part_origin, variable_origin, map_value)

        # Adding Values
        self.nearestNeighborMapper.InverseMap(variable_origin, variable_destination, True)

        if (self.GiD_output):
            self.WriteNodalResultsCustom(self.gid_io_origin, self.model_part_origin, variable_origin, output_time + 0.1)

        self.CheckValues(self.interface_sub_model_part_origin, variable_origin, map_value*2)

    def TestMapConstantVectorValues(self, output_time):
        map_value = [15.99, -2.88, 3.123]
        variable_origin = VELOCITY
        variable_destination = VELOCITY

        self.SetValuesOnNodes(self.model_part_origin, VELOCITY_X, map_value[0])
        self.SetValuesOnNodes(self.model_part_origin, VELOCITY_Y, map_value[1])
        self.SetValuesOnNodes(self.model_part_origin, VELOCITY_Z, map_value[2])

        if (self.GiD_output):
            self.WriteNodalResultsCustom(self.gid_io_origin, self.model_part_origin, variable_origin, output_time)

        # Overwriting Values
        self.nearestNeighborMapper.Map(variable_origin, variable_destination)

        if (self.GiD_output):
            self.WriteNodalResultsCustom(self.gid_io_destination, self.model_part_destination, variable_destination, output_time)

        self.CheckValues(self.interface_sub_model_part_destination, VELOCITY_X, map_value[0])
        self.CheckValues(self.interface_sub_model_part_destination, VELOCITY_Y, map_value[1])
        self.CheckValues(self.interface_sub_model_part_destination, VELOCITY_Z, map_value[2])

        # Adding Values
        self.nearestNeighborMapper.Map(variable_origin, variable_destination, True)

        if (self.GiD_output):
            self.WriteNodalResultsCustom(self.gid_io_destination, self.model_part_destination, variable_destination, output_time + 0.1)

        self.CheckValues(self.interface_sub_model_part_destination, VELOCITY_X, map_value[0]*2)
        self.CheckValues(self.interface_sub_model_part_destination, VELOCITY_Y, map_value[1]*2)
        self.CheckValues(self.interface_sub_model_part_destination, VELOCITY_Z, map_value[2]*2)

    def TestInverseMapConstantVectorValues(self, output_time):
        map_value = [1.4785, -0.88, -33.123]
        variable_origin = VELOCITY
        variable_destination = VELOCITY

        self.SetValuesOnNodes(self.model_part_destination, VELOCITY_X, map_value[0])
        self.SetValuesOnNodes(self.model_part_destination, VELOCITY_Y, map_value[1])
        self.SetValuesOnNodes(self.model_part_destination, VELOCITY_Z, map_value[2])

        if (self.GiD_output):
            self.WriteNodalResultsCustom(self.gid_io_destination, self.model_part_destination, variable_destination, output_time)

        # Overwriting Values
        self.nearestNeighborMapper.InverseMap(variable_origin, variable_destination)

        if (self.GiD_output):
            self.WriteNodalResultsCustom(self.gid_io_origin, self.model_part_origin, variable_origin, output_time)

        self.CheckValues(self.interface_sub_model_part_origin, VELOCITY_X, map_value[0])
        self.CheckValues(self.interface_sub_model_part_origin, VELOCITY_Y, map_value[1])
        self.CheckValues(self.interface_sub_model_part_origin, VELOCITY_Z, map_value[2])

        # Adding Values
        self.nearestNeighborMapper.InverseMap(variable_origin, variable_destination, True)

        if (self.GiD_output):
            self.WriteNodalResultsCustom(self.gid_io_origin, self.model_part_origin, variable_origin, output_time + 0.1)

        self.CheckValues(self.interface_sub_model_part_origin, VELOCITY_X, map_value[0]*2)
        self.CheckValues(self.interface_sub_model_part_origin, VELOCITY_Y, map_value[1]*2)
        self.CheckValues(self.interface_sub_model_part_origin, VELOCITY_Z, map_value[2]*2)


    def TestMapNonConstantScalarValues(self, output_time):
        variable_origin = PRESSURE
        variable_destination = PRESSURE

        self.SetValuesOnNodesId(self.interface_sub_model_part_origin, variable_origin)

        if (self.GiD_output):
            self.WriteNodalResultsCustom(self.gid_io_origin, self.model_part_origin, variable_origin, output_time)

        self.nearestNeighborMapper.Map(variable_origin, variable_destination)

        if (self.GiD_output):
            self.WriteNodalResultsCustom(self.gid_io_destination, self.model_part_destination, variable_destination, output_time)

        nodal_values = [12,14,13,22,25,23,42,59,60,57,76,77,75,89,93,112,113,115,114,121,124,123,148,147,154]
        self.CheckValuesId(self.interface_sub_model_part_destination, variable_destination, nodal_values)

    def TestInverseMapNonConstantScalarValues(self, output_time):
        variable_origin = PRESSURE
        variable_destination = PRESSURE

        self.SetValuesOnNodesId(self.interface_sub_model_part_destination, variable_destination)

        if (self.GiD_output):
            self.WriteNodalResultsCustom(self.gid_io_destination, self.model_part_destination, variable_destination, output_time)

        self.nearestNeighborMapper.InverseMap(variable_origin, variable_destination)

        if (self.GiD_output):
            self.WriteNodalResultsCustom(self.gid_io_origin, self.model_part_origin, variable_origin, output_time)

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


    def InitializeGiD(self):
        # Initialize GidIO
        output_file_origin = "KratosExecuteNearestNeighborMapperTest_gid_output/output_origin"
        output_file_destination = "KratosExecuteNearestNeighborMapperTest_gid_output/output_destination"

        gid_mode = GiDPostMode.GiD_PostAscii
        multifile = MultiFileFlag.MultipleFiles
        deformed_mesh_flag = WriteDeformedMeshFlag.WriteUndeformed
        write_conditions = WriteConditionsFlag.WriteConditions

        self.gid_io_origin = GidIO(output_file_origin, gid_mode, multifile, deformed_mesh_flag, write_conditions)
        self.gid_io_destination  = GidIO(output_file_destination,  gid_mode, multifile, deformed_mesh_flag, write_conditions)

        # Initialize Results Output
        self.gid_io_origin.InitializeResults(0, self.model_part_origin.GetMesh())
        self.gid_io_destination.InitializeResults( 0, self.model_part_destination.GetMesh())

        # Print original meshes
        self.write_mesh(self.model_part_origin, self.gid_io_origin)
        self.write_mesh(self.model_part_destination, self.gid_io_destination)

    def WriteNodalResultsCustom(self, gid_io, model_part, variable, output_time):
        gid_io.WriteNodalResults(variable, model_part.Nodes, output_time, 0)

    ### Function to write meshes ###
    def write_mesh(self, model_part, gid_io):
        gid_io.InitializeMesh(0)
        gid_io.WriteMesh(model_part.GetMesh())
        gid_io.FinalizeMesh()
