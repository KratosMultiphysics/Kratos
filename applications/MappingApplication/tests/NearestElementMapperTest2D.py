from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

from KratosMultiphysics import *
try: # test to import the modules for the parallel execution
    from KratosMultiphysics.mpi import *
    from KratosMultiphysics.MetisApplication import *
    from KratosMultiphysics.TrilinosApplication import *
except:
    pass
from KratosMultiphysics.MappingApplication import *

CheckForPreviousImport()

import KratosMultiphysics.KratosUnittest as KratosUnittest

class NearestElementMapperTest2D(KratosUnittest.TestCase):

    def __init__(self, gid_output):
        self.GiD_output = gid_output

        # Mdpa Input files
        input_file_origin      = "NearestElementMapperTest2D_mdpa/origin"
        input_file_destination = "NearestElementMapperTest2D_mdpa/destination"

        parameter_file = open("NearestElementMapperTest2D.json",'r')
        project_parameters = Parameters(parameter_file.read())
        project_parameters_mapper_1 = project_parameters["mapper_settings"][0]

        variable_list = [PRESSURE, TEMPERATURE, FORCE, VELOCITY]

        # test if executed in parallel
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

        self.model = Model()

        self.model_part_origin = self.partition_and_read_model_part(self.model,
                                                                     "ModelPartNameOrigin",
                                                                     input_file_origin, 3,
                                                                     variable_list,
                                                                     self.num_processors)

        self.model_part_destination = self.partition_and_read_model_part(self.model,
                                                                         "ModelPartNameDestination",
                                                                         input_file_destination, 3,
                                                                         variable_list,
                                                                         self.num_processors)

        # needed for the tests only, usually one does not need to get the submodel-parts for the mapper
        self.interface_sub_model_part_origin = self.model_part_origin.GetSubModelPart("interface_origin")
        self.interface_sub_model_part_destination = self.model_part_destination.GetSubModelPart("interface_destination")

        # Initialize Mapper
        self.mapper = MapperFactory.CreateMapper(self.model_part_origin,
                                    self.model_part_destination,
                                    project_parameters_mapper_1)

        if (self.GiD_output):
            self.InitializeGiD()

        # self.PrintValuesToPrescribe() # needed to set up the test

        self.SetPrescribedValues()


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
                                         Mapper.ADD_VALUES)

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
                                         Mapper.SWAP_SIGN)

        if (self.GiD_output):
            self.WriteNodalResultsCustom(self.gid_io_destination,
                                         self.model_part_destination,
                                         variable_destination,
                                         output_time + 0.2)

        self.CheckValues(self.interface_sub_model_part_destination,
                         variable_destination,
                         -map_value)

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
                                                Mapper.ADD_VALUES)

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
                                                Mapper.ADD_VALUES | Mapper.SWAP_SIGN)

        if (self.GiD_output):
            self.WriteNodalResultsCustom(self.gid_io_origin,
                                         self.model_part_origin,
                                         variable_origin,
                                         output_time + 0.2)

        self.CheckValues(self.interface_sub_model_part_destination,
                         variable_destination,
                         map_value)


        self.mapper.UpdateInterface(Mapper.REMESHED)

    def TestMapConstantVectorValues(self, output_time):
        map_value = Vector([15.99, -2.88, 3.123])
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
                                         Mapper.ADD_VALUES)

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
                                         Mapper.SWAP_SIGN)

        if (self.GiD_output):
            self.WriteNodalResultsCustom(self.gid_io_destination,
                                         self.model_part_destination,
                                         variable_destination,
                                         output_time + 0.2)

        self.CheckValues(self.interface_sub_model_part_destination,
                         variable_destination,
                         [-x for x in map_value])

    def TestInverseMapConstantVectorValues(self, output_time):
        map_value = Vector([1.4785, -0.88, -33.123])
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
                                                Mapper.ADD_VALUES)

        if (self.GiD_output):
            self.WriteNodalResultsCustom(self.gid_io_origin,
                                         self.model_part_origin,
                                         variable_origin,
                                         output_time + 0.1)

        self.CheckValues(self.interface_sub_model_part_origin,
                         variable_origin,
                         [2*x for x in map_value])

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

        # self.PrintMappedValues(self.interface_sub_model_part_destination, variable_destination) # needed to set up test
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

        # self.PrintMappedValues(self.interface_sub_model_part_origin, variable_origin) # needed to set up test
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

        # self.PrintMappedValues(self.interface_sub_model_part_destination, variable_destination) # needed to set up test
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

        # self.PrintMappedValues(self.interface_sub_model_part_origin, variable_origin) # needed to set up test
        self.CheckValuesPrescribed(self.interface_sub_model_part_origin,
                                   variable_origin,
                                   self.vector_values_origin_receive)

    def partition_and_read_model_part(self, current_model, model_part_name,
                                      model_part_input_file,
                                      size_domain, variable_list,
                                      number_of_partitions):
        model_part = current_model.CreateModelPart(model_part_name)
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
        value_to_prescribe = nodal_values[nodal_coords]
        if isinstance(value_to_prescribe, tuple):
            value_to_prescribe = Vector(list(value_to_prescribe))
        node.SetSolutionStepValue(variable, value_to_prescribe)


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


    def InitializeGiD(self):
        # Initialize GidIO
        if (self.parallel_execution):
            output_file_origin = "NearestElementMapperTest2D_gid_output/output_origin_r" + str(mpi.rank)
            output_file_destination = "NearestElementMapperTest2D_gid_output/output_destination_r" + str(mpi.rank)
        else:
            output_file_origin = "NearestElementMapperTest2D_gid_output/output_origin"
            output_file_destination = "NearestElementMapperTest2D_gid_output/output_destination"

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

    def PrintValuesToPrescribe(self):
        from random import uniform
        values_range = [-100, 100]
        values_precision = 4
        print("Origin ModelPart; Scalar Values")
        for node in self.interface_sub_model_part_origin.Nodes:
            nodal_coords = (node.X, node.Y, node.Z)
            value = round(uniform(values_range[0],values_range[1]),values_precision)
            print(str(nodal_coords) + " : " + str(value) + ",")

        print("\n\nDestination ModelPart; Scalar Values")
        for node in self.interface_sub_model_part_destination.Nodes:
            nodal_coords = (node.X, node.Y, node.Z)
            value = round(uniform(values_range[0],values_range[1]),values_precision)
            print(str(nodal_coords) + " : " + str(value) + ",")

        print("\n\nOrigin ModelPart; Vector Values")
        for node in self.interface_sub_model_part_origin.Nodes:
            nodal_coords = (node.X, node.Y, node.Z)
            value_x = round(uniform(values_range[0],values_range[1]),values_precision)
            value_y = round(uniform(values_range[0],values_range[1]),values_precision)
            value_z = round(uniform(values_range[0],values_range[1]),values_precision)
            values = (value_x, value_y, value_z)
            print(str(nodal_coords) + " : " + str(values) + ",")

        print("\n\nDestination ModelPart; Vector Values")
        for node in self.interface_sub_model_part_destination.Nodes:
            nodal_coords = (node.X, node.Y, node.Z)
            value_x = round(uniform(values_range[0],values_range[1]),values_precision)
            value_y = round(uniform(values_range[0],values_range[1]),values_precision)
            value_z = round(uniform(values_range[0],values_range[1]),values_precision)
            values = (value_x, value_y, value_z)
            print(str(nodal_coords) + " : " + str(values) + ",")

        err # needed to get the output

    def PrintMappedValues(self, model_part, variable):
        for node in model_part.Nodes:
            nodal_coords = (node.X, node.Y, node.Z)
            print(str(nodal_coords) + " : " + str(node.GetSolutionStepValue(variable)) + ",")
        err # needed to get the output

    def SetPrescribedValues(self):
        self.scalar_values_origin_send = {
        (0.5, 0.0, 0.0)       : 53.4779,
        (0.45652, 0.0, 0.0)   : 5.2858,
        (0.41304, 0.0, 0.0)   : 59.0475,
        (0.36957, 0.0, 0.0)   : 27.9523,
        (0.32609, 0.0, 0.0)   : -32.1865,
        (0.28261, 0.0, 0.0)   : 63.8048,
        (0.23913, 0.0, 0.0)   : 71.511,
        (0.19565, 0.0, 0.0)   : -13.1623,
        (0.15217, 0.0, 0.0)   : 95.9792,
        (0.1087, 0.0, 0.0)    : -30.2111,
        (0.06522, 0.0, 0.0)   : 30.0971,
        (0.02174, 0.0, 0.0)   : 83.4711,
        (-0.02174, 0.0, 0.0)  : -14.6857,
        (-0.06522, 0.0, 0.0)  : 38.0837,
        (-0.1087, 0.0, 0.0)   : -86.7614,
        (-0.15217, 0.0, 0.0)  : -14.0549,
        (-0.19565, 0.0, 0.0)  : 71.5862,
        (-0.23913, 0.0, 0.0)  : -27.7837,
        (-0.28261, 0.0, 0.0)  : -6.31,
        (-0.32609, 0.0, 0.0)  : 45.6769,
        (-0.36957, 0.0, 0.0)  : -14.3858,
        (-0.41304, 0.0, 0.0)  : -50.2728,
        (-0.45652, 0.0, 0.0)  : 95.5271,
        (-0.5, 0.0, 0.0)      : -94.2814 }

        self.scalar_values_destination_send = {
        (0.5, 0.0, 0.0)       : -31.0195,
        (0.42857, 0.0, 0.0)   : 4.4518,
        (0.35714, 0.0, 0.0)   : -10.4431,
        (0.28571, 0.0, 0.0)   : 45.9975,
        (0.21429, 0.0, 0.0)   : 12.1486,
        (0.14286, 0.0, 0.0)   : 58.7815,
        (0.07143, 0.0, 0.0)   : -70.8426,
        (0.0, 0.0, 0.0)       : 91.1817,
        (-0.07143, 0.0, 0.0)  : -24.2759,
        (-0.14286, 0.0, 0.0)  : -82.7191,
        (-0.21429, 0.0, 0.0)  : 11.5811,
        (-0.28571, 0.0, 0.0)  : 3.6705,
        (-0.35714, 0.0, 0.0)  : -17.4966,
        (-0.42857, 0.0, 0.0)  : -16.8597,
        (-0.5, 0.0, 0.0)      : -20.9641 }

        self.vector_values_origin_send = {
        (0.5, 0.0, 0.0)       : (29.6949, -67.9985, -70.2535),
        (0.45652, 0.0, 0.0)   : (-22.4765, -82.5025, -12.0273),
        (0.41304, 0.0, 0.0)   : (85.6208, -37.5963, -70.4871),
        (0.36957, 0.0, 0.0)   : (26.8462, 51.5921, -54.9098),
        (0.32609, 0.0, 0.0)   : (-69.2305, 96.9028, -90.0268),
        (0.28261, 0.0, 0.0)   : (52.5312, -80.2247, -15.1366),
        (0.23913, 0.0, 0.0)   : (-31.4985, -95.3868, 54.7422),
        (0.19565, 0.0, 0.0)   : (50.7086, -95.1887, -70.0518),
        (0.15217, 0.0, 0.0)   : (17.8997, -39.1578, -93.9703),
        (0.1087, 0.0, 0.0)    : (-35.6544, 91.5894, -32.0882),
        (0.06522, 0.0, 0.0)   : (-12.7211, -74.838, 89.1694),
        (0.02174, 0.0, 0.0)   : (-68.6042, -53.681, 14.9236),
        (-0.02174, 0.0, 0.0)  : (-17.2847, 48.0229, -59.9835),
        (-0.06522, 0.0, 0.0)  : (61.2927, 23.0026, 88.8573),
        (-0.1087, 0.0, 0.0)   : (85.6711, 74.239, -9.7172),
        (-0.15217, 0.0, 0.0)  : (94.5805, 32.3122, 28.2284),
        (-0.19565, 0.0, 0.0)  : (-39.6473, 1.9494, 77.8997),
        (-0.23913, 0.0, 0.0)  : (4.2095, -89.4119, 18.3684),
        (-0.28261, 0.0, 0.0)  : (77.1388, -40.5598, 96.9675),
        (-0.32609, 0.0, 0.0)  : (19.6688, 21.0438, -29.4639),
        (-0.36957, 0.0, 0.0)  : (27.1732, 46.5045, -38.9491),
        (-0.41304, 0.0, 0.0)  : (48.8117, 59.4878, 63.8178),
        (-0.45652, 0.0, 0.0)  : (-87.5449, 11.0379, 1.8602),
        (-0.5, 0.0, 0.0)      : (90.4201, -67.3753, 48.7639) }

        self.vector_values_destination_send = {
        (0.5, 0.0, 0.0)       : (-88.7441, 45.4286, 98.9649),
        (0.42857, 0.0, 0.0)   : (0.6243, 9.0312, 96.0099),
        (0.35714, 0.0, 0.0)   : (54.8269, 77.9277, -26.0114),
        (0.28571, 0.0, 0.0)   : (-49.7987, -32.3158, -50.6387),
        (0.21429, 0.0, 0.0)   : (-74.9104, 20.1111, -27.7886),
        (0.14286, 0.0, 0.0)   : (-51.1795, 47.0775, 70.399),
        (0.07143, 0.0, 0.0)   : (-2.0198, 48.6805, -83.8233),
        (0.0, 0.0, 0.0)       : (11.4188, 5.5089, -43.3743),
        (-0.07143, 0.0, 0.0)  : (-78.2683, 84.1022, 63.1146),
        (-0.14286, 0.0, 0.0)  : (35.9396, 13.3631, 4.7501),
        (-0.21429, 0.0, 0.0)  : (-23.0865, -99.0268, -90.4462),
        (-0.28571, 0.0, 0.0)  : (-46.4922, -57.032, -5.5928),
        (-0.35714, 0.0, 0.0)  : (8.463, -8.884, -13.4453),
        (-0.42857, 0.0, 0.0)  : (-38.7724, 36.9184, -1.1619),
        (-0.5, 0.0, 0.0)      : (-46.1269, 79.0725, -22.9492) }


        self.scalar_values_origin_receive = {
        (0.5, 0.0, 0.0)       : -31.0195,
        (0.45652, 0.0, 0.0)   : -9.4278,
        (0.41304, 0.0, 0.0)   : 1.2134,
        (0.36957, 0.0, 0.0)   : -7.8511,
        (0.32609, 0.0, 0.0)   : 14.0911,
        (0.28261, 0.0, 0.0)   : 44.5283,
        (0.23913, 0.0, 0.0)   : 23.9213,
        (0.19565, 0.0, 0.0)   : 24.3177,
        (0.15217, 0.0, 0.0)   : 52.7035,
        (0.1087, 0.0, 0.0)    : -3.2087,
        (0.06522, 0.0, 0.0)   : -56.7565,
        (0.02174, 0.0, 0.0)   : 41.869,
        (-0.02174, 0.0, 0.0)  : 56.0417,
        (-0.06522, 0.0, 0.0)  : -14.2382,
        (-0.1087, 0.0, 0.0)   : -54.7698,
        (-0.15217, 0.0, 0.0)  : -70.4283,
        (-0.19565, 0.0, 0.0)  : -13.027,
        (-0.23913, 0.0, 0.0)  : 8.8298,
        (-0.28261, 0.0, 0.0)  : 4.0139,
        (-0.32609, 0.0, 0.0)  : -8.2954,
        (-0.36957, 0.0, 0.0)  : -17.3858,
        (-0.41304, 0.0, 0.0)  : -16.9982,
        (-0.45652, 0.0, 0.0)  : -18.4657,
        (-0.5, 0.0, 0.0)      : -20.9641 }

        self.scalar_values_destination_receive = {
        (0.5, 0.0, 0.0)       : 53.4779,
        (0.42857, 0.0, 0.0)   : 39.8451,
        (0.35714, 0.0, 0.0)   : 10.7599,
        (0.28571, 0.0, 0.0)   : 56.9609,
        (0.21429, 0.0, 0.0)   : 23.1374,
        (0.14286, 0.0, 0.0)   : 68.9529,
        (0.07143, 0.0, 0.0)   : 21.4836,
        (0.0, 0.0, 0.0)       : 34.3927,
        (-0.07143, 0.0, 0.0)  : 20.2528,
        (-0.14286, 0.0, 0.0)  : -29.6265,
        (-0.21429, 0.0, 0.0)  : 28.986,
        (-0.28571, 0.0, 0.0)  : -2.6035,
        (-0.35714, 0.0, 0.0)  : 2.7848,
        (-0.42857, 0.0, 0.0)  : 1.8034,
        (-0.5, 0.0, 0.0)      : -94.2814 }

        self.vector_values_origin_receive = {
        (0.5, 0.0, 0.0)       : (-88.7441, 45.4286, 98.9649),
        (0.45652, 0.0, 0.0)   : (-34.3449, 23.2732, 97.1662),
        (0.41304, 0.0, 0.0)   : (12.4088, 24.0104, 69.4806),
        (0.36957, 0.0, 0.0)   : (45.3948, 65.9386, -4.7777),
        (0.32609, 0.0, 0.0)   : (9.3471, 30.0058, -36.7167),
        (0.28261, 0.0, 0.0)   : (-50.8887, -30.0402, -49.6469),
        (0.23913, 0.0, 0.0)   : (-66.1765, 1.8769, -35.7359),
        (0.19565, 0.0, 0.0)   : (-68.7177, 27.1481, -2.1661),
        (0.15217, 0.0, 0.0)   : (-54.2725, 43.5628, 57.6015),
        (0.1087, 0.0, 0.0)    : (-27.6698, 47.8441, -3.3548),
        (0.06522, 0.0, 0.0)   : (-0.8515, 44.9272, -80.3067),
        (0.02174, 0.0, 0.0)   : (7.3287, 18.6483, -55.6851),
        (-0.02174, 0.0, 0.0)  : (-15.8778, 29.4291, -10.9640),
        (-0.06522, 0.0, 0.0)  : (-70.4711, 77.2694, 53.8566),
        (-0.1087, 0.0, 0.0)   : (-18.6781, 47.1927, 32.6618),
        (-0.15217, 0.0, 0.0)  : (28.2463, -1.2855, -7.6575),
        (-0.19565, 0.0, 0.0)  : (-7.6834, -69.6981, -65.6043),
        (-0.23913, 0.0, 0.0)  : (-31.2270, -84.4209, -60.9340),
        (-0.28261, 0.0, 0.0)  : (-45.4763, -58.8548, -9.2759),
        (-0.32609, 0.0, 0.0)  : (-15.4255, -29.8135, -10.0319),
        (-0.36957, 0.0, 0.0)  : (0.2433, -0.9136, -11.3078),
        (-0.41304, 0.0, 0.0)  : (-28.5027, 26.9602, -3.8325),
        (-0.45652, 0.0, 0.0)  : (-41.6502, 53.4130, -9.6871),
        (-0.5, 0.0, 0.0)      : (-46.1269, 79.0725, -22.9492) }

        self.vector_values_destination_receive = {
        (0.5, 0.0, 0.0)       : (29.6949, -67.9985, -70.2535),
        (0.42857, 0.0, 0.0)   : (47.0111, -53.6357, -49.6067),
        (0.35714, 0.0, 0.0)   : (-0.6201, 64.5455, -64.9490),
        (0.28571, 0.0, 0.0)   : (43.8499, -67.5960, -20.4761),
        (0.21429, 0.0, 0.0)   : (15.4662, -95.2736, -16.5523),
        (0.14286, 0.0, 0.0)   : (6.4300, -11.1556, -80.7170),
        (0.07143, 0.0, 0.0)   : (-15.9965, -51.0681, 71.8509),
        (0.0, 0.0, 0.0)       : (-42.9445, -2.8291, -22.5300),
        (-0.07143, 0.0, 0.0)  : (64.7745, 30.3204, 74.7785),
        (-0.14286, 0.0, 0.0)  : (92.6724, 41.2917, 20.1016),
        (-0.21429, 0.0, 0.0)  : (-20.8458, -37.2174, 52.3785),
        (-0.28571, 0.0, 0.0)  : (73.0414, -36.1676, 87.9533),
        (-0.35714, 0.0, 0.0)  : (25.0279, 39.2258, -36.2375),
        (-0.42857, 0.0, 0.0)  : (0.1084, 42.1827, 41.6881),
        (-0.5, 0.0, 0.0)      : (90.4201, -67.3753, 48.7639) }
