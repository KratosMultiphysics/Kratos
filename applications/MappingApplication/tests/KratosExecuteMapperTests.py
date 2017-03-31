from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

from KratosMultiphysics import *
from KratosMultiphysics.MappingApplication import *

try: # test to import the modules for the parallel execution
    from KratosMultiphysics.mpi import *
    from KratosMultiphysics.MetisApplication import *
    from KratosMultiphysics.TrilinosApplication import *
except:
    pass

CheckForPreviousImport()

import KratosMultiphysics.KratosUnittest as KratosUnittest


class KratosExecuteMapperTests(KratosUnittest.TestCase):

    def __init__(self, Gid_Output):
        self.GiD_output = Gid_Output

        # Mdpa Input files
        input_file_origin      = "MapperTests_mdpa/MappingApplication_test_geometry_tri"
        input_file_destination = "MapperTests_mdpa/MappingApplication_test_geometry_quad"

        self.variable_list_scalar = [PRESSURE, TEMPERATURE]
        self.variable_list_vector = [FORCE, VELOCITY]
        
        variable_list = []
        variable_list.extend(self.variable_list_scalar)
        variable_list.extend(self.variable_list_vector)

        # check if executed in parallel
        try:
            self.num_processors = mpi.size
        except:
            self.num_processors = 1

        if (self.num_processors == 1): # serial execution
            self.parallel_execution = False
        else:
            # Partition and Read Model Parts
            variable_list.extend([PARTITION_INDEX])
            self.parallel_execution = True

        self.model_part_origin = self.partition_and_read_model_part("ModelPartNameOrigin",
                                                                     input_file_origin, 3,
                                                                     variable_list,
                                                                     self.num_processors)
        self.model_part_destination = self.partition_and_read_model_part("ModelPartNameDestination",
                                                                         input_file_destination, 3,
                                                                         variable_list,
                                                                         self.num_processors)  

    def SetUpMapper(self, FileName):
        self.ResetValuesModelParts()

        if (self.GiD_output):
            self.InitializeGiD(FileName)

        parameter_file = open("MapperTests_json/" + FileName + "_parameters.json", 'r')
        result_file    = open("MapperTests_results/" + FileName + "_results.json", 'r')
        
        project_parameters = Parameters(parameter_file.read())
        self.results =            Parameters(result_file.read())

        self.SetPrescribedValues()

        mapper_settings = project_parameters["mapper_settings"][0]

        # needed for the tests only, usually one does not need to get the submodel-parts for the mapper
        self.interface_sub_model_part_origin = self.model_part_origin.GetSubModelPart(
                                                mapper_settings["interface_submodel_part_origin"].GetString())
        self.interface_sub_model_part_destination = self.model_part_destination.GetSubModelPart(
                                                mapper_settings["interface_submodel_part_destination"].GetString())

        # Initialize Mapper
        self.mapper = MapperFactory(self.model_part_origin,
                                    self.model_part_destination,
                                    mapper_settings)
        

        # self.PrintValuesForJson() # needed to set up the test


    def TestMapConstantScalarValues(self, output_time):
        map_value = 5.123
        variable_origin = PRESSURE
        variable_destination = TEMPERATURE

        self.SetValuesOnNodes(self.model_part_origin,
                              variable_origin,
                              map_value)

        if (self.GiD_output):
            self.WriteNodalResultsCustom(self.gid_io_origin,
                                         self.model_part_origin,
                                         variable_origin,
                                         output_time)


        # Overwriting Values
        self.mapper.Map(variable_origin,
                        variable_destination)

        if (self.GiD_output):
            self.WriteNodalResultsCustom(self.gid_io_destination,
                                         self.model_part_destination,
                                         variable_destination,
                                         output_time)

        self.CheckValues(self.interface_sub_model_part_destination,
                         variable_destination, map_value)


        # Adding Values
        self.mapper.Map(variable_origin,
                        variable_destination,
                        MapperFactory.ADD_VALUES)

        if (self.GiD_output):
            self.WriteNodalResultsCustom(self.gid_io_destination,
                                         self.model_part_destination,
                                         variable_destination,
                                         output_time + 0.1)

        self.CheckValues(self.interface_sub_model_part_destination,
                         variable_destination, map_value*2)


        # Swaping the sign of the Values
        self.mapper.Map(variable_origin,
                                         variable_destination,
                                         MapperFactory.SWAP_SIGN)

        if (self.GiD_output):
            self.WriteNodalResultsCustom(self.gid_io_destination,
                                         self.model_part_destination,
                                         variable_destination,
                                         output_time + 0.2)

        self.CheckValues(self.interface_sub_model_part_destination,
                         variable_destination,
                         -map_value)


        # # Conservative Mapping
        # # Number of Nodes on Origin: 37
        # # Number of Nodes in Destination: 25
        # # => Values in Destination are multiplied with a factor of 1.48 (37/25)
        # # to conserve the sum of quantities aka conservative mapping
        # self.mapper.Map(variable_origin,
        #                                  variable_destination,
        #                                  MapperFactory.CONSERVATIVE)

        # if (self.GiD_output):
        #     self.WriteNodalResultsCustom(self.gid_io_destination,
        #                                  self.model_part_destination,
        #                                  variable_destination,
        #                                  output_time + 0.3)

        # self.CheckValues(self.interface_sub_model_part_destination,
        #                  variable_destination,
        #                  map_value*1.48)


        self.mapper.UpdateInterface()

    def TestInverseMapConstantScalarValues(self, output_time):
        map_value = -8.6647
        variable_origin = TEMPERATURE
        variable_destination = PRESSURE

        self.SetValuesOnNodes(self.model_part_destination,
                              variable_destination,
                              map_value)

        if (self.GiD_output):
            self.WriteNodalResultsCustom(self.gid_io_destination,
                                         self.model_part_destination,
                                         variable_destination,
                                         output_time)


        # Overwriting Values
        self.mapper.InverseMap(variable_origin,
                                                variable_destination)

        if (self.GiD_output):
            self.WriteNodalResultsCustom(self.gid_io_origin,
                                         self.model_part_origin,
                                         variable_origin,
                                         output_time)

        self.CheckValues(self.interface_sub_model_part_origin,
                         variable_origin,
                         map_value)


        # Adding Values
        self.mapper.InverseMap(variable_origin,
                                                variable_destination,
                                                MapperFactory.ADD_VALUES)

        if (self.GiD_output):
            self.WriteNodalResultsCustom(self.gid_io_origin,
                                         self.model_part_origin,
                                         variable_origin,
                                         output_time + 0.1)

        self.CheckValues(self.interface_sub_model_part_origin,
                         variable_origin,
                         map_value*2)


        # Swaping the sign of the Values and adding them
        self.mapper.InverseMap(variable_origin,
                                                variable_destination,
                                                MapperFactory.ADD_VALUES | MapperFactory.SWAP_SIGN)

        if (self.GiD_output):
            self.WriteNodalResultsCustom(self.gid_io_origin,
                                         self.model_part_origin,
                                         variable_origin,
                                         output_time + 0.2)

        self.CheckValues(self.interface_sub_model_part_destination,
                         variable_destination,
                         map_value)


        self.mapper.UpdateInterface(MapperFactory.REMESHED)

    def TestMapConstantVectorValues(self, output_time):
        map_value = [15.99, -2.88, 3.123]
        variable_origin = FORCE
        variable_destination = VELOCITY

        self.SetValuesOnNodes(self.model_part_origin,
                              variable_origin,
                              map_value)

        if (self.GiD_output):
            self.WriteNodalResultsCustom(self.gid_io_origin,
                                         self.model_part_origin,
                                         variable_origin,
                                         output_time)


        # Overwriting Values
        self.mapper.Map(variable_origin,
                                         variable_destination)

        if (self.GiD_output):
            self.WriteNodalResultsCustom(self.gid_io_destination,
                                         self.model_part_destination,
                                         variable_destination,
                                         output_time)

        self.CheckValues(self.interface_sub_model_part_destination,
                         variable_destination,
                         map_value)


        # Adding Values
        self.mapper.Map(variable_origin,
                                         variable_destination,
                                         MapperFactory.ADD_VALUES)

        if (self.GiD_output):
            self.WriteNodalResultsCustom(self.gid_io_destination,
                                         self.model_part_destination,
                                         variable_destination,
                                         output_time + 0.1)

        self.CheckValues(self.interface_sub_model_part_destination,
                         variable_destination,
                         [2*x for x in map_value])


        # Swaping the sign of the Values
        self.mapper.Map(variable_origin,
                                         variable_destination,
                                         MapperFactory.SWAP_SIGN)

        if (self.GiD_output):
            self.WriteNodalResultsCustom(self.gid_io_destination,
                                         self.model_part_destination,
                                         variable_destination,
                                         output_time + 0.2)

        self.CheckValues(self.interface_sub_model_part_destination,
                         variable_destination,
                         [-x for x in map_value])

    def TestInverseMapConstantVectorValues(self, output_time):
        map_value = [1.4785, -0.88, -33.123]
        variable_origin = VELOCITY
        variable_destination = FORCE

        self.SetValuesOnNodes(self.model_part_destination,
                              variable_destination,
                              map_value)

        if (self.GiD_output):
            self.WriteNodalResultsCustom(self.gid_io_destination,
                                         self.model_part_destination,
                                         variable_destination,
                                         output_time)

        # Overwriting Values
        self.mapper.InverseMap(variable_origin,
                                                variable_destination)

        if (self.GiD_output):
            self.WriteNodalResultsCustom(self.gid_io_origin,
                                         self.model_part_origin,
                                         variable_origin,
                                         output_time)

        self.CheckValues(self.interface_sub_model_part_origin,
                         variable_origin,
                         map_value)

        # Adding Values
        self.mapper.InverseMap(variable_origin,
                                                variable_destination,
                                                MapperFactory.ADD_VALUES)

        if (self.GiD_output):
            self.WriteNodalResultsCustom(self.gid_io_origin,
                                         self.model_part_origin,
                                         variable_origin,
                                         output_time + 0.1)

        self.CheckValues(self.interface_sub_model_part_origin,
                         variable_origin,
                         [2*x for x in map_value])

        # # Conservative Mapping
        # # Number of Nodes on Origin: 37
        # # Number of Nodes in Destination: 25
        # # => Values in Origin are multiplied with a factor of 0.675675676 (25/37)
        # # to conserve the sum of quantities aka conservative mapping
        # self.mapper.InverseMap(variable_origin,
        #                                         variable_destination,
        #                                         MapperFactory.CONSERVATIVE)

        # if (self.GiD_output):
        #     self.WriteNodalResultsCustom(self.gid_io_origin,
        #                                  self.model_part_origin,
        #                                  variable_origin,
        #                                  output_time + 0.2)

        # self.CheckValues(self.interface_sub_model_part_origin,
        #                  variable_origin,
        #                  [0.675675676*x for x in map_value])


    def TestMapNonConstantScalarValues(self, output_time):
        variable_origin = PRESSURE
        variable_destination = TEMPERATURE

        self.SetValuesOnNodesPrescribed(self.interface_sub_model_part_origin,
                                        variable_origin,
                                        self.scalar_values_origin_send)

        if (self.GiD_output):
            self.WriteNodalResultsCustom(self.gid_io_origin,
                                         self.model_part_origin,
                                         variable_origin,
                                         output_time)

        self.mapper.Map(variable_origin,
                                         variable_destination)

        if (self.GiD_output):
            self.WriteNodalResultsCustom(self.gid_io_destination,
                                         self.model_part_destination,
                                         variable_destination,
                                         output_time)

        # self.PrintMappedValues(self.interface_sub_model_part_destination, variable_destination,
        #                        "Destination_receive Scalar") # needed to set up test

        self.CheckValuesPrescribed(self.interface_sub_model_part_destination,
                                   variable_destination,
                                   self.scalar_values_destination_receive)

    def TestInverseMapNonConstantScalarValues(self, output_time):
        variable_origin = TEMPERATURE
        variable_destination = PRESSURE

        self.SetValuesOnNodesPrescribed(self.interface_sub_model_part_destination,
                                        variable_destination,
                                        self.scalar_values_destination_send)

        if (self.GiD_output):
            self.WriteNodalResultsCustom(self.gid_io_destination,
                                         self.model_part_destination,
                                         variable_destination,
                                         output_time)

        self.mapper.InverseMap(variable_origin,
                                                variable_destination)

        if (self.GiD_output):
            self.WriteNodalResultsCustom(self.gid_io_origin,
                                         self.model_part_origin,
                                         variable_origin,
                                         output_time)

        self.PrintMappedValues(self.interface_sub_model_part_origin, variable_origin,
                               "Origin_receive Scalar") # needed to set up test

        self.CheckValuesPrescribed(self.interface_sub_model_part_origin,
                                   variable_origin,
                                   self.scalar_values_origin_receive)

    def TestMapNonConstantVectorValues(self, output_time):
        variable_origin = FORCE
        variable_destination = VELOCITY

        self.SetValuesOnNodesPrescribed(self.interface_sub_model_part_origin,
                                        variable_origin,
                                        self.vector_values_origin_send)

        if (self.GiD_output):
            self.WriteNodalResultsCustom(self.gid_io_origin,
                                         self.model_part_origin,
                                         variable_origin,
                                         output_time)

        self.mapper.Map(variable_origin,
                                         variable_destination)

        if (self.GiD_output):
            self.WriteNodalResultsCustom(self.gid_io_destination,
                                         self.model_part_destination,
                                         variable_destination,
                                         output_time)

        self.PrintMappedValues(self.interface_sub_model_part_destination, variable_destination,
                               "Destination_receive Vector") # needed to set up test

        self.CheckValuesPrescribed(self.interface_sub_model_part_destination,
                                   variable_destination,
                                   self.vector_values_destination_receive)

    def TestInverseMapNonConstantVectorValues(self, output_time):
        variable_origin = VELOCITY
        variable_destination = FORCE

        self.SetValuesOnNodesPrescribed(self.interface_sub_model_part_destination,
                                        variable_destination,
                                        self.vector_values_destination_send)

        if (self.GiD_output):
            self.WriteNodalResultsCustom(self.gid_io_destination,
                                         self.model_part_destination,
                                         variable_destination,
                                         output_time)

        self.mapper.InverseMap(variable_origin,
                                                variable_destination)

        if (self.GiD_output):
            self.WriteNodalResultsCustom(self.gid_io_origin,
                                         self.model_part_origin,
                                         variable_origin,
                                         output_time)

        self.PrintMappedValues(self.interface_sub_model_part_origin, variable_origin,
                               "Origin_receive Vector") # needed to set up test

        self.CheckValuesPrescribed(self.interface_sub_model_part_origin,
                                   variable_origin,
                                   self.vector_values_origin_receive)

    def partition_and_read_model_part(self, model_part_name,
                                      model_part_input_file,
                                      size_domain, variable_list,
                                      number_of_partitions):
        model_part = ModelPart(model_part_name)
        for variable in variable_list:
            model_part.AddNodalSolutionStepVariable(variable)

        if (number_of_partitions > 1):
            if (mpi.size > 1):
                if (mpi.rank == 0):
                    model_part_io = ReorderConsecutiveModelPartIO(model_part_input_file)

                    partitioner = MetisDivideHeterogeneousInputProcess(
                        model_part_io,
                        number_of_partitions,
                        size_domain,
                        0, # verbosity, set to 1 for more detailed output
                        True)

                    partitioner.Execute()

                mpi.world.barrier()
                model_part_input_file = model_part_input_file + "_" + str(mpi.rank)

        model_part_io = ModelPartIO(model_part_input_file)
        model_part_io.ReadModelPart(model_part)

        if (number_of_partitions > 1):
            MPICommSetup = SetMPICommunicatorProcess(model_part)
            MPICommSetup.Execute()

            ParallelFillComm = ParallelFillCommunicator(model_part.GetRootModelPart())
            ParallelFillComm.Execute()

        model_part.ProcessInfo.SetValue(DOMAIN_SIZE, size_domain)
        model_part.SetBufferSize(1)

        return model_part

    def ResetValuesModelParts(self):
        for node in self.model_part_origin.Nodes:
            for variable in self.variable_list_scalar:
                node.SetSolutionStepValue(variable, 0.0)
            for variable in self.variable_list_vector:
                node.SetSolutionStepValue(variable, (0.0, 0.0, 0.0))

        for node in self.model_part_destination.Nodes:
            for variable in self.variable_list_scalar:
                node.SetSolutionStepValue(variable, 0.0)
            for variable in self.variable_list_vector:
                node.SetSolutionStepValue(variable, (0.0, 0.0, 0.0))

    def SetValuesOnNodes(self, model_part, variable, value):
        for node in model_part.Nodes:
            if (self.parallel_execution):
                if (node.GetSolutionStepValue(PARTITION_INDEX) == mpi.rank):
                    self.SetValuesOnNodesExec(node, variable, value)
            else:
                self.SetValuesOnNodesExec(node, variable, value)

    def SetValuesOnNodesExec(self, node, variable, value):
        node.SetSolutionStepValue(variable, value)

    def SetValuesOnNodesPrescribed(self, model_part, variable, nodal_values):
        for node in model_part.Nodes:
            if (self.parallel_execution):
                if (node.GetSolutionStepValue(PARTITION_INDEX) == mpi.rank):
                    self.SetValuesOnNodesPrescribedExec(node, variable, nodal_values)
            else:
                self.SetValuesOnNodesPrescribedExec(node, variable, nodal_values)

    def SetValuesOnNodesPrescribedExec(self, node, variable, nodal_values):
        nodal_coords = (node.X, node.Y, node.Z)
        node.SetSolutionStepValue(variable, nodal_values[nodal_coords])


    def CheckValues(self, model_part, variable, value_mapped):
        for node in model_part.Nodes:
            if (self.parallel_execution):
                if (node.GetSolutionStepValue(PARTITION_INDEX) == mpi.rank):
                    self.CheckValuesExec(node, variable, value_mapped)
            else:
                self.CheckValuesExec(node, variable, value_mapped)

    def CheckValuesExec(self, node, variable, value_mapped):
        value_expected = node.GetSolutionStepValue(variable)
        if (isinstance(value_mapped, float) or isinstance(value_mapped, int)): # Variable is a scalar
            self.assertAlmostEqual(value_mapped,value_expected,4)
        else: # Variable is a vector
            self.assertAlmostEqualVector(value_mapped,value_expected)


    def CheckValuesPrescribed(self, model_part, variable, nodal_values):
        for node in model_part.Nodes:
            if (self.parallel_execution):
                if (node.GetSolutionStepValue(PARTITION_INDEX) == mpi.rank):
                    self.CheckValuesPrescribedExec(node, variable, nodal_values)
            else:
                self.CheckValuesPrescribedExec(node, variable, nodal_values)

    def CheckValuesPrescribedExec(self, node, variable, nodal_values):
        value_mapped = node.GetSolutionStepValue(variable)
        nodal_coords = (node.X, node.Y, node.Z)
        value_expected = nodal_values[nodal_coords]
        if (isinstance(value_mapped, float) or isinstance(value_mapped, int)): # Variable is a scalar
            self.assertAlmostEqual(value_mapped,value_expected,4)


    def assertAlmostEqualVector(self, values_mapped, values_expected):
        for i in range(0,3):
            self.assertAlmostEqual(values_mapped[i],values_expected[i],4)


    def InitializeGiD(self, FileName):
        # Initialize GidIO
        output_file_origin = "MapperTests_gid_output/output_" + FileName + "_origin"
        output_file_destination = "MapperTests_gid_output/output_" + FileName + "_destination"

        if (self.parallel_execution):
            output_file_origin += "_r" + str(mpi.rank)
            output_file_destination += "_r" + str(mpi.rank)

        gid_mode = GiDPostMode.GiD_PostAscii
        multifile = MultiFileFlag.MultipleFiles
        deformed_mesh_flag = WriteDeformedMeshFlag.WriteUndeformed
        write_conditions = WriteConditionsFlag.WriteConditions

        self.gid_io_origin = GidIO(output_file_origin, gid_mode, multifile,
                                   deformed_mesh_flag, write_conditions)
        self.gid_io_destination = GidIO(output_file_destination, gid_mode, multifile,
                                        deformed_mesh_flag, write_conditions)

        # Initialize Results Output
        self.gid_io_origin.InitializeResults(0, self.model_part_origin.GetMesh())
        self.gid_io_destination.InitializeResults( 0, self.model_part_destination.GetMesh())

        # Print original meshes
        self.write_mesh(self.model_part_origin, self.gid_io_origin)
        self.write_mesh(self.model_part_destination, self.gid_io_destination)

    def WriteNodalResultsCustom(self, gid_io, model_part, variable, output_time):
        if (self.parallel_execution):
            gid_io.WriteNodalResults(variable, model_part.GetCommunicator().LocalMesh().Nodes, output_time, 0)
        else:
            gid_io.WriteNodalResults(variable, model_part.Nodes, output_time, 0)

    ### Function to write meshes ###
    def write_mesh(self, model_part, gid_io):
        gid_io.InitializeMesh(0)
        gid_io.WriteMesh(model_part.GetMesh())
        gid_io.FinalizeMesh()

    def FinalizeGiD(self):
        # Finalize GidIO
        self.gid_io_origin.FinalizeResults()
        self.gid_io_destination.FinalizeResults()

    def PrintValuesForJson(self):
        coordinates_origin_x = []
        coordinates_origin_y = []
        coordinates_origin_z = []
        for node in self.interface_sub_model_part_origin.Nodes:
            coordinates_origin_x.append(node.X)
            coordinates_origin_y.append(node.Y)
            coordinates_origin_z.append(node.Z)

        coordinates_destination_x = []
        coordinates_destination_y = []
        coordinates_destination_z = []
        for node in self.interface_sub_model_part_destination.Nodes:
            coordinates_destination_x.append(node.X)
            coordinates_destination_y.append(node.Y)
            coordinates_destination_z.append(node.Z)

        num_nodes_origin = len(coordinates_origin_x)
        print(num_nodes_origin)
        num_nodes_destination = len(coordinates_destination_x)
        print(num_nodes_destination)

        print("\n\nCopy into JSON:\n")
        print("{")
        print("    \"Origin_Coordinates\": {")
        print("        \"X\":", coordinates_origin_x, ",")
        print("        \"Y\":", coordinates_origin_y, ",")
        print("        \"Z\":", coordinates_origin_z)
        print("    },")
        print("    \"Destination_Coordinates\": {")
        print("        \"X\":", coordinates_destination_x, ",")
        print("        \"Y\":", coordinates_destination_y, ",")
        print("        \"Z\":", coordinates_destination_z)
        print("    },")
        print("    \"Origin_send\": {")
        print("        \"Scalar\":", self.CreateRandomValues(num_nodes_origin), ",")
        print("        \"Vector_X\":", self.CreateRandomValues(num_nodes_origin), ",")
        print("        \"Vector_Y\":", self.CreateRandomValues(num_nodes_origin), ",")
        print("        \"Vector_Z\":", self.CreateRandomValues(num_nodes_origin))
        print("    },")
        print("    \"Destination_send\": {")
        print("        \"Scalar\":", self.CreateRandomValues(num_nodes_destination), ",")
        print("        \"Vector_X\":", self.CreateRandomValues(num_nodes_destination), ",")
        print("        \"Vector_Y\":", self.CreateRandomValues(num_nodes_destination), ",")
        print("        \"Vector_Z\":", self.CreateRandomValues(num_nodes_destination))
        print("    },")
        print("    \"Origin_receive\": {")
        print("        \"Scalar\": [],")
        print("        \"Vector_X\": [],")
        print("        \"Vector_Y\": [],")
        print("        \"Vector_Z\": []")
        print("    },")
        print("    \"Destination_receive\": {")
        print("        \"Scalar\": [],")
        print("        \"Vector_X\": [],")
        print("        \"Vector_Y\": [],")
        print("        \"Vector_Z\": []")
        print("    }")
        print("}")

        err # needed to get the output

    def CreateRandomValues(self, NumValues):
        from random import uniform
        values_range = [-100, 100]
        values_precision = 4
        values = []
        for i in range(NumValues) :
            values.append(round(uniform(values_range[0],values_range[1]),values_precision))

        return values

    def PrintMappedValues(self, model_part, variable, InfoString):
        if ("Vector" in InfoString): # Vector Variable
            values_x = []
            values_y = []
            values_z = []
            for node in model_part.Nodes:
                value = node.GetSolutionStepValue(variable)
                values_x.append(value[0])
                values_y.append(value[1])
                values_z.append(value[2])
            print("\n\n Values for:", InfoString, "X, Y, Z")
            print(values_x)
            print(values_y)
            print(values_z)

        else: # Scalar Variable
            values = []
            for node in model_part.Nodes:
                values.append(node.GetSolutionStepValue(variable))
            print("\n\n Values for:", InfoString)
            print(values)
        err # needed to get the output

    def SetPrescribedValues(self):
        self.scalar_values_origin_send = {}
        self.scalar_values_destination_send = {}
        self.scalar_values_origin_receive = {}
        self.scalar_values_destination_receive = {}

        self.vector_values_origin_send = {}
        self.vector_values_destination_send = {}
        self.vector_values_origin_receive = {}
        self.vector_values_destination_receive = {}

        # values_origin_send_scalar = [ Results["Origin_send"]["Scalar"] ]
        # values_destination_send_scalar = [ Results["Destination_send"]["Scalar"] ]
        # values_origin_receive_scalar = [ Results["Origin_receive"]["Scalar"] ]
        # values_destination_receive_scalar = [ Results["Destination_receive"]["Scalar"] ]

        # values_origin_send_vector = [ Results["Origin_send"]["Vector_X"], 
        #                               Results["Origin_send"]["Vector_Y"], 
        #                               Results["Origin_send"]["Vector_Z"] ]
        
        # values_destination_send_vector = [ Results["Destination_send"]["Vector_X"], 
        #                                    Results["Destination_send"]["Vector_Y"], 
        #                                    Results["Destination_send"]["Vector_Z"] ]

        # values_origin_receive_vector = [ Results["Origin_receive"]["Vector_X"], 
        #                                  Results["Origin_receive"]["Vector_Y"], 
        #                                  Results["Origin_receive"]["Vector_Z"] ]

        # values_destination_revceive_vector = [ Results["Destination_receive"]["Vector_X"], 
        #                                        Results["Destination_receive"]["Vector_Y"], 
        #                                        Results["Destination_receive"]["Vector_Z"] ]

        self.AssignPrescribedValues("Origin",
                                    "send", 
                                    "Scalar",
                                    self.scalar_values_origin_send)

        self.AssignPrescribedValues("Origin",
                                    "send", 
                                    "Vector",
                                    self.vector_values_origin_send)

        self.AssignPrescribedValues("Destination",
                                    "send", 
                                    "Scalar",
                                    self.scalar_values_destination_send)

        self.AssignPrescribedValues("Destination",
                                    "send", 
                                    "Vector",
                                    self.vector_values_destination_send)

        # self.AssignPrescribedScalarValues(self.coordinates_destination,
        #                                   values_destination_send_scalar,
        #                                   self.scalar_values_destination_send)
                                          
        # # self.AssignPrescribedScalarValues(self.coordinates_origin, 
        # #                                   values_origin_receive_scalar,
        # #                                   self.scalar_values_origin_receive)

        # # self.AssignPrescribedScalarValues(self.coordinates_destination, 
        # #                                   values_destination_receive_scalar,
        # #                                   self.scalar_values_destination_receive)


        # self.AssignPrescribedScalarValues(self.coordinates_origin, 
        #                                   values_origin_send_vector,
        #                                   self.vector_values_origin_send)

        # self.AssignPrescribedScalarValues(self.coordinates_destination, 
        #                                   values_destination_send_vector,
        #                                   self.vector_values_destination_send)
                                          
        # self.AssignPrescribedScalarValues(self.coordinates_origin, 
        #                                   values_origin_receive_vector,
        #                                   self.vector_values_origin_receive)

        # self.AssignPrescribedScalarValues(self.coordinates_destination, 
        #                                   values_destination_receive_vector,
        #                                   self.vector_values_destination_receive)
        # err

    def AssignPrescribedValues(self, side_interface, direction, value_type, dictionary):
        coordinates_x = self.results[side_interface + "_Coordinates"]["X"]
        coordinates_y = self.results[side_interface + "_Coordinates"]["Y"]
        coordinates_z = self.results[side_interface + "_Coordinates"]["Z"]

        values = 0
        if (value_type == "Scalar"):
            values = self.results[side_interface + "_" + direction]["Scalar"]
        elif(value_type == "Vector"):
            values = [ self.results[side_interface + "_" + direction]["Vector_X"],
                       self.results[side_interface + "_" + direction]["Vector_Y"],
                       self.results[side_interface + "_" + direction]["Vector_Z"] ]
        else:
            raise("ERROR: wrong value_type")

        i = 0
        for x in coordinates_x:
            # Retreive Coordinates
            coords = (coordinates_x[i].GetDouble(), 
                      coordinates_y[i].GetDouble(), 
                      coordinates_z[i].GetDouble())

            # Retreive Values
            if (value_type == "Scalar"):
                value = values[i].GetDouble()
            elif(value_type == "Vector"):
                value = (values[0][i].GetDouble(),
                         values[1][i].GetDouble(),
                         values[2][i].GetDouble())
            
            dictionary.update({coords : value})

            i += 1



    # def SetModelPartCoordinates(self):
    #     # Origin ModelPart (Triangles / Tetrahedrals)
    #     self.coordinates_origin = [
    #     (0.5, 0.5, 0.3), (0.5, 0.20106, 0.3), (0.5, 0.5, 0.0), (0.16709, 0.5, 0.3), 
    #     (0.19321, 0.29439, 0.14288), (0.5, 0.20106, 0.0), (0.16709, 0.5, 0.0), 
    #     (0.25444, 0.11511, 0.3), (0.03593, 0.25692, 0.3), (0.5, -0.05492, 0.3), 
    #     (0.26229, 0.07307, 0.0), (0.1209, 0.04321, 0.21031), (0.5, -0.05492, 0.0), 
    #     (-0.16581, 0.5, 0.3), (0.00041, 0.15728, 0.0), (-0.16581, 0.5, 0.0), 
    #     (0.01629, -0.06596, 0.3), (0.25366, -0.21527, 0.3), (0.1, -0.10072, 0.0), 
    #     (0.5, -0.28585, 0.3), (-0.22958, 0.11452, 0.3), (0.5, -0.28585, 0.0), 
    #     (-0.31205, 0.359, 0.1), (-0.0001, -0.20352, 0.15055), (0.25, -0.29263, 0.0), 
    #     (-0.1, -0.10072, 0.0), (-0.26112, 0.07206, 0.0), (0.24884, -0.39894, 0.12901), 
    #     (0.5, -0.5, 0.3), (-0.5, 0.5, 0.3), (0.0, -0.31427, 0.0), (-0.24285, -0.21189, 0.3), 
    #     (0.25699, -0.5, 0.3), (-0.00169, -0.38579, 0.10309), (-0.5, 0.20068, 0.3), 
    #     (0.5, -0.5, 0.0), (-0.5, 0.5, 0.0), (-0.33063, -0.13518, 0.16015), (0.3, -0.5, 0.0), 
    #     (-0.5, 0.20068, 0.0), (0.1, -0.5, 0.0), (-0.00577, -0.5, 0.3), (-0.25, -0.29315, 0.0), 
    #     (-0.5, -0.05656, 0.3), (-0.22592, -0.39278, 0.12864), (-0.5, -0.05656, 0.0), 
    #     (-0.1, -0.5, 0.0), (-0.37761, -0.37663, 0.1513), (-0.25983, -0.5, 0.3), 
    #     (-0.5, -0.28735, 0.3), (-0.5, -0.28735, 0.0), (-0.3, -0.5, 0.0), (-0.5, -0.5, 0.3), 
    #     (-0.5, -0.5, 0.0)]

    #     # Destination ModelPart (Quads / Hexas)
    #     self.coordinates_destination = [
    #     (0.4, 0.4, 0.0), (0.4, 0.4, 0.15), (0.175, 0.4, 0.0), (0.175, 0.4, 0.15), 
    #     (0.4, 0.4, 0.3), (0.4, 0.1, 0.0), (0.4, 0.1, 0.15), (0.175, 0.1, 0.0), 
    #     (0.175, 0.4, 0.3), (0.175, 0.1, 0.15), (0.4, 0.1, 0.3), (-0.05, 0.4, 0.0), 
    #     (-0.05, 0.4, 0.15), (0.175, 0.1, 0.3), (-0.05, 0.4, 0.3), (-0.05, 0.1, 0.0), 
    #     (-0.05, 0.1, 0.15), (0.4, -0.2, 0.0), (0.4, -0.2, 0.15), (-0.05, 0.1, 0.3), 
    #     (0.175, -0.2, 0.0), (0.175, -0.2, 0.15), (0.4, -0.2, 0.3), (-0.275, 0.4, 0.0), 
    #     (-0.275, 0.4, 0.15), (0.175, -0.2, 0.3), (-0.275, 0.4, 0.3), (-0.275, 0.1, 0.0), 
    #     (-0.05, -0.2, 0.0), (-0.275, 0.1, 0.15), (-0.05, -0.2, 0.15), (-0.275, 0.1, 0.3), 
    #     (-0.05, -0.2, 0.3), (0.4, -0.5, 0.0), (-0.5, 0.4, 0.0), (-0.275, -0.2, 0.0), 
    #     (-0.5, 0.4, 0.15), (0.4, -0.5, 0.15), (-0.275, -0.2, 0.15), (0.175, -0.5, 0.0), 
    #     (0.175, -0.5, 0.15), (0.4, -0.5, 0.3), (-0.5, 0.4, 0.3), (-0.5, 0.1, 0.0), 
    #     (-0.275, -0.2, 0.3), (-0.5, 0.1, 0.15), (0.175, -0.5, 0.3), (-0.5, 0.1, 0.3), 
    #     (-0.05, -0.5, 0.0), (-0.05, -0.5, 0.15), (-0.05, -0.5, 0.3), (-0.5, -0.2, 0.0), 
    #     (-0.5, -0.2, 0.15), (-0.5, -0.2, 0.3), (-0.275, -0.5, 0.0), (-0.275, -0.5, 0.15), 
    #     (-0.275, -0.5, 0.3), (-0.5, -0.5, 0.0), (-0.5, -0.5, 0.15), (-0.5, -0.5, 0.3)]


      
  
  
  
    # # def SetPrescribedValues(self):
    # #     self.scalar_values_origin_send = {
    # #     (-0.5, -0.5, 0.0) : 74.8612,
    # #     (-0.3002, -0.5, 0.0) : 71.2017,
    # #     (-0.5, -0.364, 0.0) : -80.2386,
    # #     (-0.32243, -0.21657, 0.0) : 8.5714,
    # #     (-0.199, -0.5, 0.0) : -81.6182,
    # #     (-0.15359, -0.33, 0.0) : -85.2777,
    # #     (-0.5, -0.17, 0.0) : -94.3847,
    # #     (-0.15359, -0.1004, 0.0) : 70.1815,
    # #     (-0.32679, -0.00077, 0.0) : -95.3669,
    # #     (0.01401, -0.35733, 0.0) : 41.267,
    # #     (0.173, -0.5, 0.0) : 77.2532,
    # #     (-0.5, 0.188, 0.0) : -58.9351,
    # #     (-0.15359, 0.1888, 0.0) : -75.2223,
    # #     (0.18145, -0.31476, 0.0) : -1.8728,
    # #     (0.01262, -0.04, 0.0) : -88.0393,
    # #     (-0.32293, 0.23667, 0.0) : -0.6181,
    # #     (0.19212, -0.177, 0.0) : 99.6869,
    # #     (-0.5, 0.355, 0.0) : -7.3811,
    # #     (0.3568, -0.5, 0.0) : 29.8251,
    # #     (-0.15559, 0.355, 0.0) : -69.2206,
    # #     (0.33541, -0.21888, 0.0) : 57.3482,
    # #     (0.03398, 0.21551, 0.0) : 61.6527,
    # #     (0.19244, 0.1066, 0.0) : 88.1457,
    # #     (0.32613, 0.012, 0.0) : 92.5446,
    # #     (0.5, -0.5, 0.0) : 46.5348,
    # #     (-0.5, 0.5, 0.0) : -65.896,
    # #     (-0.31004, 0.5, 0.0) : -72.4578,
    # #     (0.5, -0.325, 0.0) : 65.5213,
    # #     (0.19885, 0.29804, 0.0) : -79.2115,
    # #     (0.34366, 0.15961, 0.0) : 0.3687,
    # #     (-0.10041, 0.5, 0.0) : 0.5273,
    # #     (0.5, -0.10561, 0.0) : 17.4257,
    # #     (0.152, 0.5, 0.0) : 67.2235,
    # #     (0.5, 0.19874, 0.0) : -91.6093,
    # #     (0.3321, 0.5, 0.0) : -70.3876,
    # #     (0.5, 0.30456, 0.0) : 87.3811,
    # #     (0.5, 0.5, 0.0) : 47.4853 }

    # #     self.scalar_values_destination_send = {
    # #     (-0.5, -0.5, 0.0) : -29.2379,
    # #     (-0.5, -0.246, 0.0) : -98.8008,
    # #     (-0.262, -0.5, 0.0) : -21.5213,
    # #     (-0.28857, -0.14832, 0.0) : 78.5838,
    # #     (-0.5, 0.06104, 0.0) : 86.4897,
    # #     (0.0553, -0.5, 0.0) : 33.4003,
    # #     (-0.04361, -0.28681, 0.0) : 60.6458,
    # #     (-0.08274, -0.01745, 0.0) : 44.3008,
    # #     (-0.27185, 0.1518, 0.0) : -14.2951,
    # #     (0.15858, -0.32561, 0.0) : 45.8605,
    # #     (-0.5, 0.2592, 0.0) : -38.8276,
    # #     (0.2511, -0.5, 0.0) : 20.6577,
    # #     (0.19743, -0.1191, 0.0) : -70.3052,
    # #     (-0.14373, 0.25658, 0.0) : 11.9008,
    # #     (0.19981, 0.11831, 0.0) : -10.9708,
    # #     (0.5, -0.5, 0.0) : -70.6218,
    # #     (-0.5, 0.5, 0.0) : -71.0203,
    # #     (0.5, -0.2425, 0.0) : 57.4228,
    # #     (-0.25129, 0.5, 0.0) : -99.7596,
    # #     (0.11949, 0.31894, 0.0) : 21.5272,
    # #     (0.5, 0.0052, 0.0) : -22.0285,
    # #     (0.0124, 0.5, 0.0) : 35.9431,
    # #     (0.5, 0.249, 0.0) : 87.495,
    # #     (0.25186, 0.5, 0.0) : 93.018,
    # #     (0.5, 0.5, 0.0) : 13.4472 }

    # #     self.vector_values_origin_send = {
    # #     (-0.5, -0.5, 0.0) : (-82.6987, 95.4827, 89.3428),
    # #     (-0.3002, -0.5, 0.0) : (60.7885, 89.6246, -15.217),
    # #     (-0.5, -0.364, 0.0) : (60.5731, -87.3619, 14.2162),
    # #     (-0.32243, -0.21657, 0.0) : (94.3739, -47.977, 84.8945),
    # #     (-0.199, -0.5, 0.0) : (-45.466, 10.6292, 13.0954),
    # #     (-0.15359, -0.33, 0.0) : (70.8957, -79.9669, 10.832),
    # #     (-0.5, -0.17, 0.0) : (97.4695, 81.8736, -88.5944),
    # #     (-0.15359, -0.1004, 0.0) : (-94.7962, 70.3255, 40.8216),
    # #     (-0.32679, -0.00077, 0.0) : (-86.3453, 30.0437, -63.3355),
    # #     (0.01401, -0.35733, 0.0) : (81.9431, 44.044, -76.4111),
    # #     (0.173, -0.5, 0.0) : (-80.3916, -32.4842, 67.057),
    # #     (-0.5, 0.188, 0.0) : (-46.661, 23.6675, 89.1151),
    # #     (-0.15359, 0.1888, 0.0) : (47.3517, 34.2461, -56.1094),
    # #     (0.18145, -0.31476, 0.0) : (68.6203, 49.2641, 2.576),
    # #     (0.01262, -0.04, 0.0) : (-95.6574, 80.6006, -20.3202),
    # #     (-0.32293, 0.23667, 0.0) : (35.495, -89.1297, 72.2261),
    # #     (0.19212, -0.177, 0.0) : (70.5812, -40.4829, 97.703),
    # #     (-0.5, 0.355, 0.0) : (-26.6771, -45.2197, -43.927),
    # #     (0.3568, -0.5, 0.0) : (-33.4444, -59.0702, -23.3442),
    # #     (-0.15559, 0.355, 0.0) : (-62.3037, -60.0014, 11.5194),
    # #     (0.33541, -0.21888, 0.0) : (-36.7822, -74.8982, 94.9569),
    # #     (0.03398, 0.21551, 0.0) : (87.1369, 59.4913, -78.8707),
    # #     (0.19244, 0.1066, 0.0) : (-31.8417, -14.0344, -9.4779),
    # #     (0.32613, 0.012, 0.0) : (-58.5903, 63.2413, 80.5193),
    # #     (0.5, -0.5, 0.0) : (-43.5106, -57.763, 98.5749),
    # #     (-0.5, 0.5, 0.0) : (-30.8904, 50.4188, 89.1966),
    # #     (-0.31004, 0.5, 0.0) : (93.0028, -53.8615, -34.9477),
    # #     (0.5, -0.325, 0.0) : (0.9265, -40.9413, -36.5583),
    # #     (0.19885, 0.29804, 0.0) : (70.0601, 63.2991, -11.9064),
    # #     (0.34366, 0.15961, 0.0) : (-28.9694, -26.4793, -6.908),
    # #     (-0.10041, 0.5, 0.0) : (83.3125, -14.6739, 23.429),
    # #     (0.5, -0.10561, 0.0) : (63.7904, -21.4233, 83.3817),
    # #     (0.152, 0.5, 0.0) : (77.2986, -28.2953, 35.1405),
    # #     (0.5, 0.19874, 0.0) : (20.9458, -28.2064, -7.7992),
    # #     (0.3321, 0.5, 0.0) : (12.9683, 66.6925, 28.39),
    # #     (0.5, 0.30456, 0.0) : (78.6603, -82.0786, -0.9975),
    # #     (0.5, 0.5, 0.0) : (17.363, -14.019, 98.6192) }

    # #     self.vector_values_destination_send = {
    # #     (-0.5, -0.5, 0.0) : (-56.4463, -25.9128, 65.2663),
    # #     (-0.5, -0.246, 0.0) : (15.5472, -50.5599, -6.8589),
    # #     (-0.262, -0.5, 0.0) : (-42.618, -77.049, 15.075),
    # #     (-0.28857, -0.14832, 0.0) : (-45.3971, 12.1528, -76.4542),
    # #     (-0.5, 0.06104, 0.0) : (76.7786, -77.7307, 82.6615),
    # #     (0.0553, -0.5, 0.0) : (-75.3713, 59.8855, -26.5626),
    # #     (-0.04361, -0.28681, 0.0) : (80.9893, -54.946, 9.0543),
    # #     (-0.08274, -0.01745, 0.0) : (72.7381, -47.1258, 11.7604),
    # #     (-0.27185, 0.1518, 0.0) : (52.1867, -89.2781, -25.5366),
    # #     (0.15858, -0.32561, 0.0) : (6.3482, -88.3025, -67.2723),
    # #     (-0.5, 0.2592, 0.0) : (73.4525, 56.9959, -57.0295),
    # #     (0.2511, -0.5, 0.0) : (-30.4067, 83.2038, 82.9062),
    # #     (0.19743, -0.1191, 0.0) : (-95.7395, 6.1749, 29.4679),
    # #     (-0.14373, 0.25658, 0.0) : (-44.646, -11.6643, 36.3388),
    # #     (0.19981, 0.11831, 0.0) : (8.9531, 48.6582, 22.0265),
    # #     (0.5, -0.5, 0.0) : (-63.2004, -89.0804, -47.276),
    # #     (-0.5, 0.5, 0.0) : (-62.2399, 78.8062, -43.6925),
    # #     (0.5, -0.2425, 0.0) : (-61.2615, -96.1305, 62.5688),
    # #     (-0.25129, 0.5, 0.0) : (-52.5028, 7.931, 66.3195),
    # #     (0.11949, 0.31894, 0.0) : (-43.234, -97.6892, -38.1162),
    # #     (0.5, 0.0052, 0.0) : (-9.6424, 3.8129, -3.9237),
    # #     (0.0124, 0.5, 0.0) : (99.2372, 38.2263, 92.0563),
    # #     (0.5, 0.249, 0.0) : (-6.5423, 16.3678, 0.9326),
    # #     (0.25186, 0.5, 0.0) : (-30.2543, 19.2205, 10.2341),
    # #     (0.5, 0.5, 0.0) : (91.3727, -50.0157, -27.6178) }


    # #     self.scalar_values_origin_receive = {
    # #     (-0.5, -0.5, 0.0) : -29.2379,
    # #     (-0.3002, -0.5, 0.0) : -21.5213,
    # #     (-0.5, -0.364, 0.0) : -98.8008,
    # #     (-0.32243, -0.21657, 0.0) : 78.5838,
    # #     (-0.199, -0.5, 0.0) : -21.5213,
    # #     (-0.15359, -0.33, 0.0) : 60.6458,
    # #     (-0.5, -0.17, 0.0) : -98.8008,
    # #     (-0.15359, -0.1004, 0.0) : 44.3008,
    # #     (-0.32679, -0.00077, 0.0) : 78.5838,
    # #     (0.01401, -0.35733, 0.0) : 60.6458,
    # #     (0.173, -0.5, 0.0) : 20.6577,
    # #     (-0.5, 0.188, 0.0) : -38.8276,
    # #     (-0.15359, 0.1888, 0.0) : 11.9008,
    # #     (0.18145, -0.31476, 0.0) : 45.8605,
    # #     (0.01262, -0.04, 0.0) : 44.3008,
    # #     (-0.32293, 0.23667, 0.0) : -14.2951,
    # #     (0.19212, -0.177, 0.0) : -70.3052,
    # #     (-0.5, 0.355, 0.0) : -38.8276,
    # #     (0.3568, -0.5, 0.0) : 20.6577,
    # #     (-0.15559, 0.355, 0.0) : 11.9008,
    # #     (0.33541, -0.21888, 0.0) : 57.4228,
    # #     (0.03398, 0.21551, 0.0) : 21.5272,
    # #     (0.19244, 0.1066, 0.0) : -10.9708,
    # #     (0.32613, 0.012, 0.0) : -10.9708,
    # #     (0.5, -0.5, 0.0) : -70.6218,
    # #     (-0.5, 0.5, 0.0) : -71.0203,
    # #     (-0.31004, 0.5, 0.0) : -99.7596,
    # #     (0.5, -0.325, 0.0) : 57.4228,
    # #     (0.19885, 0.29804, 0.0) : 21.5272,
    # #     (0.34366, 0.15961, 0.0) : -10.9708,
    # #     (-0.10041, 0.5, 0.0) : 35.9431,
    # #     (0.5, -0.10561, 0.0) : -22.0285,
    # #     (0.152, 0.5, 0.0) : 93.018,
    # #     (0.5, 0.19874, 0.0) : 87.495,
    # #     (0.3321, 0.5, 0.0) : 93.018,
    # #     (0.5, 0.30456, 0.0) : 87.495,
    # #     (0.5, 0.5, 0.0) : 13.4472 }

    # #     self.scalar_values_destination_receive = {
    # #     (-0.5, -0.5, 0.0) : 74.8612,
    # #     (-0.5, -0.246, 0.0) : -94.3847,
    # #     (-0.262, -0.5, 0.0) : 71.2017,
    # #     (-0.28857, -0.14832, 0.0) : 8.5714,
    # #     (-0.5, 0.06104, 0.0) : -58.9351,
    # #     (0.0553, -0.5, 0.0) : 77.2532,
    # #     (-0.04361, -0.28681, 0.0) : 41.267,
    # #     (-0.08274, -0.01745, 0.0) : -88.0393,
    # #     (-0.27185, 0.1518, 0.0) : -0.6181,
    # #     (0.15858, -0.32561, 0.0) : -1.8728,
    # #     (-0.5, 0.2592, 0.0) : -58.9351,
    # #     (0.2511, -0.5, 0.0) : 77.2532,
    # #     (0.19743, -0.1191, 0.0) : 99.6869,
    # #     (-0.14373, 0.25658, 0.0) : -75.2223,
    # #     (0.19981, 0.11831, 0.0) : 88.1457,
    # #     (0.5, -0.5, 0.0) : 46.5348,
    # #     (-0.5, 0.5, 0.0) : -65.896,
    # #     (0.5, -0.2425, 0.0) : 65.5213,
    # #     (-0.25129, 0.5, 0.0) : -72.4578,
    # #     (0.11949, 0.31894, 0.0) : -79.2115,
    # #     (0.5, 0.0052, 0.0) : 17.4257,
    # #     (0.0124, 0.5, 0.0) : 0.5273,
    # #     (0.5, 0.249, 0.0) : -91.6093,
    # #     (0.25186, 0.5, 0.0) : -70.3876,
    # #     (0.5, 0.5, 0.0) : 47.4853 }

    # #     self.vector_values_origin_receive = {
    # #     (-0.5, -0.5, 0.0) : (-56.4463,-25.9128,65.2663),
    # #     (-0.3002, -0.5, 0.0) : (-42.618,-77.049,15.075),
    # #     (-0.5, -0.364, 0.0) : (15.5472,-50.5599,-6.8589),
    # #     (-0.32243, -0.21657, 0.0) : (-45.3971,12.1528,-76.4542),
    # #     (-0.199, -0.5, 0.0) : (-42.618,-77.049,15.075),
    # #     (-0.15359, -0.33, 0.0) : (80.9893,-54.946,9.0543),
    # #     (-0.5, -0.17, 0.0) : (15.5472,-50.5599,-6.8589),
    # #     (-0.15359, -0.1004, 0.0) : (72.7381,-47.1258,11.7604),
    # #     (-0.32679, -0.00077, 0.0) : (-45.3971,12.1528,-76.4542),
    # #     (0.01401, -0.35733, 0.0) : (80.9893,-54.946,9.0543),
    # #     (0.173, -0.5, 0.0) : (-30.4067,83.2038,82.9062),
    # #     (-0.5, 0.188, 0.0) : (73.4525,56.9959,-57.0295),
    # #     (-0.15359, 0.1888, 0.0) : (-44.646,-11.6643,36.3388),
    # #     (0.18145, -0.31476, 0.0) : (6.3482,-88.3025,-67.2723),
    # #     (0.01262, -0.04, 0.0) : (72.7381,-47.1258,11.7604),
    # #     (-0.32293, 0.23667, 0.0) : (52.1867,-89.2781,-25.5366),
    # #     (0.19212, -0.177, 0.0) : (-95.7395,6.1749,29.4679),
    # #     (-0.5, 0.355, 0.0) : (73.4525,56.9959,-57.0295),
    # #     (0.3568, -0.5, 0.0) : (-30.4067,83.2038,82.9062),
    # #     (-0.15559, 0.355, 0.0) : (-44.646,-11.6643,36.3388),
    # #     (0.33541, -0.21888, 0.0) : (-61.2615,-96.1305,62.5688),
    # #     (0.03398, 0.21551, 0.0) : (-43.234,-97.6892,-38.1162),
    # #     (0.19244, 0.1066, 0.0) : (8.9531,48.6582,22.0265),
    # #     (0.32613, 0.012, 0.0) : (8.9531,48.6582,22.0265),
    # #     (0.5, -0.5, 0.0) : (-63.2004,-89.0804,-47.276),
    # #     (-0.5, 0.5, 0.0) : (-62.2399,78.8062,-43.6925),
    # #     (-0.31004, 0.5, 0.0) : (-52.5028,7.931,66.3195),
    # #     (0.5, -0.325, 0.0) : (-61.2615,-96.1305,62.5688),
    # #     (0.19885, 0.29804, 0.0) : (-43.234,-97.6892,-38.1162),
    # #     (0.34366, 0.15961, 0.0) : (8.9531,48.6582,22.0265),
    # #     (-0.10041, 0.5, 0.0) : (99.2372,38.2263,92.0563),
    # #     (0.5, -0.10561, 0.0) : (-9.6424,3.8129,-3.9237),
    # #     (0.152, 0.5, 0.0) : (-30.2543,19.2205,10.2341),
    # #     (0.5, 0.19874, 0.0) : (-6.5423,16.3678,0.9326),
    # #     (0.3321, 0.5, 0.0) : (-30.2543,19.2205,10.2341),
    # #     (0.5, 0.30456, 0.0) : (-6.5423,16.3678,0.9326),
    # #     (0.5, 0.5, 0.0) : (91.3727,-50.0157,-27.6178) }

    # #     self.vector_values_destination_receive = {
    # #     (-0.5, -0.5, 0.0) : (-82.6987,95.4827,89.3428),
    # #     (-0.5, -0.246, 0.0) : (97.4695,81.8736,-88.5944),
    # #     (-0.262, -0.5, 0.0) : (60.7885,89.6246,-15.217),
    # #     (-0.28857, -0.14832, 0.0) : (94.3739,-47.977,84.8945),
    # #     (-0.5, 0.06104, 0.0) : (-46.661,23.6675,89.1151),
    # #     (0.0553, -0.5, 0.0) : (-80.3916,-32.4842,67.057),
    # #     (-0.04361, -0.28681, 0.0) : (81.9431,44.044,-76.4111),
    # #     (-0.08274, -0.01745, 0.0) : (-95.6574,80.6006,-20.3202),
    # #     (-0.27185, 0.1518, 0.0) : (35.495,-89.1297,72.2261),
    # #     (0.15858, -0.32561, 0.0) : (68.6203,49.2641,2.576),
    # #     (-0.5, 0.2592, 0.0) : (-46.661,23.6675,89.1151),
    # #     (0.2511, -0.5, 0.0) : (-80.3916,-32.4842,67.057),
    # #     (0.19743, -0.1191, 0.0) : (70.5812,-40.4829,97.703),
    # #     (-0.14373, 0.25658, 0.0) : (47.3517,34.2461,-56.1094),
    # #     (0.19981, 0.11831, 0.0) : (-31.8417,-14.0344,-9.4779),
    # #     (0.5, -0.5, 0.0) : (-43.5106,-57.763,98.5749),
    # #     (-0.5, 0.5, 0.0) : (-30.8904,50.4188,89.1966),
    # #     (0.5, -0.2425, 0.0) : (0.9265,-40.9413,-36.5583),
    # #     (-0.25129, 0.5, 0.0) : (93.0028,-53.8615,-34.9477),
    # #     (0.11949, 0.31894, 0.0) : (70.0601,63.2991,-11.9064),
    # #     (0.5, 0.0052, 0.0) : (63.7904,-21.4233,83.3817),
    # #     (0.0124, 0.5, 0.0) : (83.3125,-14.6739,23.429),
    # #     (0.5, 0.249, 0.0) : (20.9458,-28.2064,-7.7992),
    # #     (0.25186, 0.5, 0.0) : (12.9683,66.6925,28.39),
    # #     (0.5, 0.5, 0.0) : (17.363,-14.019,98.6192) }
