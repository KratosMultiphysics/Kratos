from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

import KratosMultiphysics

class MapperTestsBase(object):
    '''This is the baseclass for the Mapper-Tests
    It prepares the ModelParts and other commonly used things
    The main purpose is to only do the reading & partitioning of the
    ModelParts once for all instances
    '''
    # called only once for this class, opposed of setUp()
    @classmethod
    def setUpClass(cls):
        # Mdpa Input files
        cls.input_file_origin      = "MapperTests_mdpa/MappingApplication_test_geometry_tri"
        cls.input_file_destination = "MapperTests_mdpa/MappingApplication_test_geometry_quad"

        cls.model_part_origin = KratosMultiphysics.ModelPart("origin")
        cls.model_part_destination = KratosMultiphysics.ModelPart("destination")

        # list of variables involved in the Mapper-Tests
        cls.variables_list_scalar = [
            KratosMultiphysics.PRESSURE,
            KratosMultiphysics.TEMPERATURE
        ]
        cls.variables_list_vector = [
            KratosMultiphysics.FORCE,
            KratosMultiphysics.VELOCITY
        ]

        for model_part in (cls.model_part_origin, cls.model_part_destination):
            model_part.ProcessInfo.SetValue(
                KratosMultiphysics.DOMAIN_SIZE, 3) # needed for the partitioner!
            for variable in cls.variables_list_scalar:
                model_part.AddNodalSolutionStepVariable(variable)
            for variable in cls.variables_list_vector:
                model_part.AddNodalSolutionStepVariable(variable)

        cls._ImportModelPart()

    @classmethod
    def _ImportModelPart(cls):
        '''In an MPI-test this function performs the partitioning and the reading
        of the ModelParts
        In the base-class only reads the ModelParts
        '''
        model_part_io = KratosMultiphysics.ModelPartIO(
            cls.input_file_origin).ReadModelPart(
            cls.model_part_origin)
        model_part_io = KratosMultiphysics.ModelPartIO(
            cls.input_file_destination).ReadModelPart(
            cls.model_part_destination)

    def setUp(self):
        '''This function resets the nodal values in the ModelParts that might exist
        from previous tests, both historical and non-historical
        '''
        zero_vector = KratosMultiphysics.Vector([0.0, 0.0, 0.0])
        for model_part in (self.model_part_origin, self.model_part_destination):
            for variable in self.variables_list_scalar:
                for node in model_part.GetCommunicator().LocalMesh().Nodes:
                    node.SetSolutionStepValue(variable, 0.0)
                    node.SetValue(variable, 0.0)
            for variable in self.variables_list_vector:
                for node in model_part.GetCommunicator().LocalMesh().Nodes:
                    node.SetSolutionStepValue(variable, zero_vector)
                    node.SetValue(variable, zero_vector)

class BaseMapperTests(MapperTestsBase):
    '''This is the baseclass for the Mapper-Tests
    It defines methods that are to be tested both in serial and in MPI
    '''
    def test_NearestNeighborMapper_line(self):
        mapper_settings = KratosMultiphysics.Parameters("""{
            "mapper_type": "nearest_neighbor",
            "interface_submodel_part_origin": "LineLoad3D_mapping_line_tri",
            "interface_submodel_part_destination": "LineLoad3D_mapping_line_quad"
        }""")
        values_file_name = "nearest_neighbor_line"

        self.__ExecuteMapperTests(mapper_settings, values_file_name)

    def test_NearestNeighborMapper_surface(self):
        mapper_settings = KratosMultiphysics.Parameters("""{
            "mapper_type": "nearest_neighbor",
            "interface_submodel_part_origin": "SurfaceLoad3D_mapping_surface_tri",
            "interface_submodel_part_destination": "SurfaceLoad3D_mapping_surface_quad"
        }""")
        values_file_name = "nearest_neighbor_surface"

        self.__ExecuteMapperTests(mapper_settings, values_file_name)

    def test_NearestNeighborMapper_volume(self):
        mapper_settings = KratosMultiphysics.Parameters("""{
            "mapper_type": "nearest_neighbor"
        }""")
        values_file_name = "nearest_volume_volume"

        self.__ExecuteMapperTests(mapper_settings, values_file_name)

    def test_NearestElementMapper_line(self):
        mapper_settings = KratosMultiphysics.Parameters("""{
            "mapper_type": "nearest_element",
            "interface_submodel_part_origin": "LineLoad3D_mapping_line_tri",
            "interface_submodel_part_destination": "LineLoad3D_mapping_line_quad"
        }""")
        values_file_name = "nearest_element_line"

        self.__ExecuteMapperTests(mapper_settings, values_file_name)

    def test_NearestElementMapper_surface(self):
        mapper_settings = KratosMultiphysics.Parameters("""{
            "mapper_type": "nearest_element",
            "interface_submodel_part_origin": "SurfaceLoad3D_mapping_surface_tri",
            "interface_submodel_part_destination": "SurfaceLoad3D_mapping_surface_quad"
        }""")
        values_file_name = "nearest_element_surface"

        self.__ExecuteMapperTests(mapper_settings, values_file_name)

    def _test_NearestElementMapper_volume(self): # TODO Implement
        mapper_settings = KratosMultiphysics.Parameters("""{
            "mapper_type": "nearest_element"
        }""")
        values_file_name = "nearest_element_volume"

        self.__ExecuteMapperTests(mapper_settings, values_file_name)

    def __ExecuteMapperTests(self, mapper_settings, values_file_name):
        # Saving the ModelParts for which the values are checked
        if mapper_settings.Has("interface_submodel_part_origin"):
            self.interface_model_part_origin = self.model_part_origin.GetSubModelPart(
                mapper_settings["interface_submodel_part_origin"].GetString())
        else:
            self.interface_model_part_origin = self.model_part_origin
        if mapper_settings.Has("interface_submodel_part_destination"):
            self.interface_model_part_destination = self.model_part_destination.GetSubModelPart(
                mapper_settings["interface_submodel_part_destination"].GetString())
        else:
            self.interface_model_part_destination = self.model_part_destination

        self.mapper = self._CreateMapper(mapper_settings)
        self.__ReadValuesFiles(values_file_name)

        self.__MapConstantScalarValues()
        self.__InverseMapConstantScalarValues()

        self.__MapConstantVectorValues()
        self.__InverseMapConstantVectorValues()

        self.__MapNonConstantScalarValues()
        self.__InverseMapNonConstantScalarValues()

        self.__MapNonConstantVectorValues()
        self.__InverseMapNonConstantVectorValues()

        self.__MapConservative()

    def __ReadValuesFiles(self, values_file_name):
        pass

    def __SetUniformValuesOnNodes(self, model_part, variable, value):
        for node in model_part.GetCommunicator().LocalMesh().Nodes:
            node.SetSolutionStepValue(variable, value)

    def __CheckUniformValues(self, model_part, variable, value):
        if type(variable) == KratosMultiphysics.Array1DVariable3:
            for node in model_part.GetCommunicator().LocalMesh().Nodes:
                value_mapped = node.GetSolutionStepValue(variable)
                for i in range(0,3):
                    self.assertAlmostEqual(value_mapped[i],value[i])
        else:
            for node in model_part.GetCommunicator().LocalMesh().Nodes:
                value_mapped = node.GetSolutionStepValue(variable)
                self.assertAlmostEqual(value_mapped,value)

    def __MapConstantScalarValues(self):
        map_value = 5.123521
        variable_origin      = self.variables_list_scalar[0]
        variable_destination = self.variables_list_scalar[1]

        self.__SetUniformValuesOnNodes(
            self.model_part_origin,
            variable_origin,
            map_value)

        # Overwriting Values
        self.mapper.Map(variable_origin,
                        variable_destination)

        self.__CheckUniformValues(
            self.interface_model_part_destination,
            variable_destination,
            map_value)

    def __InverseMapConstantScalarValues(self):
        map_value = -1145.12352
        variable_origin      = self.variables_list_scalar[1]
        variable_destination = self.variables_list_scalar[1]

        self.__SetUniformValuesOnNodes(
            self.model_part_destination,
            variable_destination,
            map_value)

        # Overwriting Values
        self.mapper.InverseMap(variable_origin,
                               variable_destination)

        self.__CheckUniformValues(
            self.interface_model_part_origin,
            variable_origin,
            map_value)

    def __MapConstantVectorValues(self):
        map_value = KratosMultiphysics.Vector([1478.0445, -300250.77801, 12580.123065])
        variable_origin      = self.variables_list_vector[0]
        variable_destination = self.variables_list_vector[1]

        self.__SetUniformValuesOnNodes(
            self.model_part_origin,
            variable_origin,
            map_value)

        # Overwriting Values
        self.mapper.Map(variable_origin,
                        variable_destination)

        self.__CheckUniformValues(
            self.interface_model_part_destination,
            variable_destination,
            map_value)

    def __InverseMapConstantVectorValues(self):
        map_value = KratosMultiphysics.Vector([147448.099445, -300.778801, 120.123065])
        variable_origin      = self.variables_list_vector[1]
        variable_destination = self.variables_list_vector[1]

        self.__SetUniformValuesOnNodes(
            self.model_part_destination,
            variable_destination,
            map_value)

        # Overwriting Values
        self.mapper.InverseMap(variable_origin,
                               variable_destination)

        self.__CheckUniformValues(
            self.interface_model_part_origin,
            variable_origin,
            map_value)

    def __MapNonConstantScalarValues(self):
        pass

    def __InverseMapNonConstantScalarValues(self):
        pass

    def __MapNonConstantVectorValues(self):
        pass

    def __InverseMapNonConstantVectorValues(self):
        pass

    def __MapConservative(self):
        '''This function check if conservative mapping works properly,
        i.e. mapping with the transpose of the mapping matrix
        '''
        pass



# from KratosMultiphysics import *
# try: # test to import the modules for the parallel execution
#     from KratosMultiphysics.mpi import *
#     from KratosMultiphysics.MetisApplication import *
#     from KratosMultiphysics.TrilinosApplication import *
# except:
#     pass
# from KratosMultiphysics.MappingApplication import *

# CheckForPreviousImport()



# class KratosExecuteMapperTests(KratosUnittest.TestCase):

#     # called only once for this class, opposed of setUp()
#     @classmethod
#     def tearDownClass(cls):
#         # try:
#         err
#         mdpa_extension = "_%s_mdpa" % mpi.rank
#         tri_file_name = "MapperTests_mdpa/MappingApplication_test_geometry_tri_"
#         quad_file_name = "MapperTests_mdpa/MappingApplication_test_geometry_quad_"
#         kratos_utils.DeleteFileIfExisting(tri_file_name  + mdpa_extension)
#         kratos_utils.DeleteFileIfExisting(quad_file_name + mdpa_extension)
#         # except:
#             # pass

#     # called only once for this class, opposed of setUp()
#     def tearDown(self):
#         # try:
#         err
#         mdpa_extension = "_%s_mdpa" % mpi.rank
#         tri_file_name = "MapperTests_mdpa/MappingApplication_test_geometry_tri_"
#         quad_file_name = "MapperTests_mdpa/MappingApplication_test_geometry_quad_"
#         kratos_utils.DeleteFileIfExisting(tri_file_name  + mdpa_extension)
#         kratos_utils.DeleteFileIfExisting(quad_file_name + mdpa_extension)
#         # except:
#             # pass

#     def __init__(self, GidOutput, set_up_test_1, set_up_test_2):
#         self.GiD_output = GidOutput
#         self.set_up_test_1 = set_up_test_1
#         self.set_up_test_2 = set_up_test_2

#         # Mdpa Input files
#         input_file_origin      = "MapperTests_mdpa/MappingApplication_test_geometry_tri"
#         input_file_destination = "MapperTests_mdpa/MappingApplication_test_geometry_quad"

#         self.variable_list_scalar = [PRESSURE, TEMPERATURE]
#         self.variable_list_vector = [FORCE, VELOCITY]

#         variable_list = []
#         variable_list.extend(self.variable_list_scalar)
#         variable_list.extend(self.variable_list_vector)

#         # check if executed in parallel
#         try:
#             num_processors = mpi.size
#         except:
#             num_processors = 1

#         if (num_processors == 1): # serial execution
#             self.parallel_execution = False
#         else:
#             # Partition and Read Model Parts
#             variable_list.extend([PARTITION_INDEX])
#             self.parallel_execution = True

#         self.model_part_origin = self.partition_and_read_model_part("ModelPartNameOrigin",
#                                                                      input_file_origin, 3,
#                                                                      variable_list,
#                                                                      num_processors)
#         self.model_part_destination = self.partition_and_read_model_part("ModelPartNameDestination",
#                                                                          input_file_destination, 3,
#                                                                          variable_list,
#                                                                          num_processors)

#     def SetUpMapper(self, file_name):
#         self.ResetValuesModelParts()

#         if (self.GiD_output):
#             self.InitializeGiD(file_name)

#         parameter_file_name = "MapperTests_json/"    + file_name + "_parameters.json"
#         result_file_name    = "MapperTests_results/" + file_name + "_results.json"

#         try: # to read project parameters file
#             parameter_file = open(parameter_file_name, 'r')
#             project_parameters = Parameters(parameter_file.read())
#             mapper_settings = project_parameters["mapper_settings"][0]
#         except:
#             raise("Project Parameter JSON File \"", parameter_file_name, "\" could not be read")

#         results_read = False
#         try: # to read the result ifle
#             result_file  = open(result_file_name, 'r')
#             self.results = Parameters(result_file.read())
#             results_read = True
#         except:
#             print("Warning: Result JSON File \"", result_file_name, "\" could not be read")

#         # needed for the tests only, usually one does not need to get the submodel-parts for the mapper
#         self.interface_sub_model_part_origin = self.model_part_origin.GetSubModelPart(
#                                                 mapper_settings["interface_submodel_part_origin"].GetString())
#         self.interface_sub_model_part_destination = self.model_part_destination.GetSubModelPart(
#                                                 mapper_settings["interface_submodel_part_destination"].GetString())

#         # Initialize Mapper
#         if self.parallel_execution:
#             fct_ptr = MapperFactory.CreateMPIMapper
#             print("Creating an MPI Mapper")
#         else:
#             fct_ptr = MapperFactory.CreateMapper
#             print("Creating a serial Mapper")
#         self.mapper = fct_ptr(self.model_part_origin,
#                               self.model_part_destination,
#                               mapper_settings)

#         if (self.set_up_test_1):
#             self.PrintValuesForJson() # needed to set up the test

#         if (results_read):
#             self.SetPrescribedValues()
#         else:
#             raise("Result JSON File \"", result_file_name, "\" could not be read")

#     ##### Testing Functions
#     def TestMapConstantScalarValues(self, output_time):
#         map_value = 5.123
#         variable_origin = PRESSURE
#         variable_destination = TEMPERATURE

#         self.SetValuesOnNodes(self.model_part_origin,
#                               variable_origin,
#                               map_value)

#         if (self.GiD_output):
#             self.WriteNodalResultsCustom(self.gid_io_origin,
#                                          self.model_part_origin,
#                                          variable_origin,
#                                          output_time)


#         # Overwriting Values
#         self.mapper.Map(variable_origin,
#                         variable_destination)

#         if (self.GiD_output):
#             self.WriteNodalResultsCustom(self.gid_io_destination,
#                                          self.model_part_destination,
#                                          variable_destination,
#                                          output_time)

#         self.CheckValues(self.interface_sub_model_part_destination,
#                          variable_destination, map_value)


#         # Adding Values
#         self.mapper.Map(variable_origin,
#                         variable_destination,
#                         Mapper.ADD_VALUES)

#         if (self.GiD_output):
#             self.WriteNodalResultsCustom(self.gid_io_destination,
#                                          self.model_part_destination,
#                                          variable_destination,
#                                          output_time + 0.1)

#         self.CheckValues(self.interface_sub_model_part_destination,
#                          variable_destination, map_value*2)


#         # Swaping the sign of the Values
#         self.mapper.Map(variable_origin,
#                         variable_destination,
#                         Mapper.SWAP_SIGN)

#         if (self.GiD_output):
#             self.WriteNodalResultsCustom(self.gid_io_destination,
#                                          self.model_part_destination,
#                                          variable_destination,
#                                          output_time + 0.2)

#         self.CheckValues(self.interface_sub_model_part_destination,
#                          variable_destination,
#                          -map_value)


#         # Swaping the sign of the Values
#         self.mapper.Map(variable_origin,
#                         variable_destination,
#                         Mapper.ADD_VALUES | Mapper.SWAP_SIGN)

#         if (self.GiD_output):
#             self.WriteNodalResultsCustom(self.gid_io_destination,
#                                          self.model_part_destination,
#                                          variable_destination,
#                                          output_time + 0.3)

#         self.CheckValues(self.interface_sub_model_part_destination,
#                          variable_destination,
#                          -map_value*2)


#         # # Conservative Mapping
#         # # Number of Nodes on Origin: 37
#         # # Number of Nodes in Destination: 25
#         # # => Values in Destination are multiplied with a factor of 1.48 (37/25)
#         # # to conserve the sum of quantities aka conservative mapping
#         # self.mapper.Map(variable_origin,
#         #                                  variable_destination,
#         #                                  Mapper.CONSERVATIVE)

#         # if (self.GiD_output):
#         #     self.WriteNodalResultsCustom(self.gid_io_destination,
#         #                                  self.model_part_destination,
#         #                                  variable_destination,
#         #                                  output_time + 0.3)

#         # self.CheckValues(self.interface_sub_model_part_destination,
#         #                  variable_destination,
#         #                  map_value*1.48)


#         self.mapper.UpdateInterface()

#     def TestInverseMapConstantScalarValues(self, output_time):
#         map_value = -8.6647
#         variable_origin = TEMPERATURE
#         variable_destination = PRESSURE

#         self.SetValuesOnNodes(self.model_part_destination,
#                               variable_destination,
#                               map_value)

#         if (self.GiD_output):
#             self.WriteNodalResultsCustom(self.gid_io_destination,
#                                          self.model_part_destination,
#                                          variable_destination,
#                                          output_time)


#         # Overwriting Values
#         self.mapper.InverseMap(variable_origin,
#                                variable_destination)

#         if (self.GiD_output):
#             self.WriteNodalResultsCustom(self.gid_io_origin,
#                                          self.model_part_origin,
#                                          variable_origin,
#                                          output_time)

#         self.CheckValues(self.interface_sub_model_part_origin,
#                          variable_origin,
#                          map_value)


#         # Adding Values
#         self.mapper.InverseMap(variable_origin,
#                                variable_destination,
#                                Mapper.ADD_VALUES)

#         if (self.GiD_output):
#             self.WriteNodalResultsCustom(self.gid_io_origin,
#                                          self.model_part_origin,
#                                          variable_origin,
#                                          output_time + 0.1)

#         self.CheckValues(self.interface_sub_model_part_origin,
#                          variable_origin,
#                          map_value*2)


#         # Swaping the sign of the Values and adding them
#         self.mapper.InverseMap(variable_origin,
#                                variable_destination,
#                                Mapper.ADD_VALUES | Mapper.SWAP_SIGN)

#         if (self.GiD_output):
#             self.WriteNodalResultsCustom(self.gid_io_origin,
#                                          self.model_part_origin,
#                                          variable_origin,
#                                          output_time + 0.2)

#         self.CheckValues(self.interface_sub_model_part_destination,
#                          variable_destination,
#                          map_value)


#         self.mapper.UpdateInterface(Mapper.REMESHED)

#     def TestMapConstantVectorValues(self, output_time):
#         map_value = Vector([15.99, -2.88, 3.123])
#         variable_origin = FORCE
#         variable_destination = VELOCITY

#         self.SetValuesOnNodes(self.model_part_origin,
#                               variable_origin,
#                               map_value)

#         if (self.GiD_output):
#             self.WriteNodalResultsCustom(self.gid_io_origin,
#                                          self.model_part_origin,
#                                          variable_origin,
#                                          output_time)


#         # Overwriting Values
#         self.mapper.Map(variable_origin,
#                         variable_destination)

#         if (self.GiD_output):
#             self.WriteNodalResultsCustom(self.gid_io_destination,
#                                          self.model_part_destination,
#                                          variable_destination,
#                                          output_time)

#         self.CheckValues(self.interface_sub_model_part_destination,
#                          variable_destination,
#                          map_value)


#         # Adding Values
#         self.mapper.Map(variable_origin,
#                         variable_destination,
#                         Mapper.ADD_VALUES)

#         if (self.GiD_output):
#             self.WriteNodalResultsCustom(self.gid_io_destination,
#                                          self.model_part_destination,
#                                          variable_destination,
#                                          output_time + 0.1)

#         self.CheckValues(self.interface_sub_model_part_destination,
#                          variable_destination,
#                          [2*x for x in map_value])


#         # Swaping the sign of the Values
#         self.mapper.Map(variable_origin,
#                         variable_destination,
#                         Mapper.SWAP_SIGN)

#         if (self.GiD_output):
#             self.WriteNodalResultsCustom(self.gid_io_destination,
#                                          self.model_part_destination,
#                                          variable_destination,
#                                          output_time + 0.2)

#         self.CheckValues(self.interface_sub_model_part_destination,
#                          variable_destination,
#                          [-x for x in map_value])

#         self.mapper.UpdateInterface(0.05)

#     def TestInverseMapConstantVectorValues(self, output_time):
#         map_value = Vector([1.4785, -0.88, -33.123])
#         variable_origin = VELOCITY
#         variable_destination = FORCE

#         self.SetValuesOnNodes(self.model_part_destination,
#                               variable_destination,
#                               map_value)

#         if (self.GiD_output):
#             self.WriteNodalResultsCustom(self.gid_io_destination,
#                                          self.model_part_destination,
#                                          variable_destination,
#                                          output_time)

#         # Overwriting Values
#         self.mapper.InverseMap(variable_origin,
#                                variable_destination)

#         if (self.GiD_output):
#             self.WriteNodalResultsCustom(self.gid_io_origin,
#                                          self.model_part_origin,
#                                          variable_origin,
#                                          output_time)

#         self.CheckValues(self.interface_sub_model_part_origin,
#                          variable_origin,
#                          map_value)

#         # Adding Values
#         self.mapper.InverseMap(variable_origin,
#                                variable_destination,
#                                Mapper.ADD_VALUES)

#         if (self.GiD_output):
#             self.WriteNodalResultsCustom(self.gid_io_origin,
#                                          self.model_part_origin,
#                                          variable_origin,
#                                          output_time + 0.1)

#         self.CheckValues(self.interface_sub_model_part_origin,
#                          variable_origin,
#                          [2*x for x in map_value])

#         # # Conservative Mapping
#         # # Number of Nodes on Origin: 37
#         # # Number of Nodes in Destination: 25
#         # # => Values in Origin are multiplied with a factor of 0.675675676 (25/37)
#         # # to conserve the sum of quantities aka conservative mapping
#         # self.mapper.InverseMap(variable_origin,
#         #                                         variable_destination,
#         #                                         Mapper.CONSERVATIVE)

#         # if (self.GiD_output):
#         #     self.WriteNodalResultsCustom(self.gid_io_origin,
#         #                                  self.model_part_origin,
#         #                                  variable_origin,
#         #                                  output_time + 0.2)

#         # self.CheckValues(self.interface_sub_model_part_origin,
#         #                  variable_origin,
#         #                  [0.675675676*x for x in map_value])


#     def TestMapNonConstantScalarValues(self, output_time):
#         variable_origin = PRESSURE
#         variable_destination = TEMPERATURE

#         self.SetValuesOnNodesPrescribed(self.interface_sub_model_part_origin,
#                                         variable_origin,
#                                         self.scalar_values_origin_send)

#         if (self.GiD_output):
#             self.WriteNodalResultsCustom(self.gid_io_origin,
#                                          self.model_part_origin,
#                                          variable_origin,
#                                          output_time)

#         self.mapper.Map(variable_origin,
#                         variable_destination)

#         if (self.GiD_output):
#             self.WriteNodalResultsCustom(self.gid_io_destination,
#                                          self.model_part_destination,
#                                          variable_destination,
#                                          output_time)

#         if (self.set_up_test_2):
#             self.PrintMappedValues(self.interface_sub_model_part_destination,
#                                    variable_destination,
#                                    "Destination_receive Scalar")
#         else:
#             self.CheckValuesPrescribed(self.interface_sub_model_part_destination,
#                                        variable_destination,
#                                        self.scalar_values_destination_receive)

#     def TestInverseMapNonConstantScalarValues(self, output_time):
#         variable_origin = TEMPERATURE
#         variable_destination = PRESSURE

#         self.SetValuesOnNodesPrescribed(self.interface_sub_model_part_destination,
#                                         variable_destination,
#                                         self.scalar_values_destination_send)

#         if (self.GiD_output):
#             self.WriteNodalResultsCustom(self.gid_io_destination,
#                                          self.model_part_destination,
#                                          variable_destination,
#                                          output_time)

#         self.mapper.InverseMap(variable_origin,
#                                variable_destination)

#         if (self.GiD_output):
#             self.WriteNodalResultsCustom(self.gid_io_origin,
#                                          self.model_part_origin,
#                                          variable_origin,
#                                          output_time)

#         if (self.set_up_test_2):
#             self.PrintMappedValues(self.interface_sub_model_part_origin,
#                                    variable_origin,
#                                    "Origin_receive Scalar")
#         else:
#             self.CheckValuesPrescribed(self.interface_sub_model_part_origin,
#                                        variable_origin,
#                                        self.scalar_values_origin_receive)

#     def TestMapNonConstantVectorValues(self, output_time):
#         variable_origin = FORCE
#         variable_destination = VELOCITY

#         self.SetValuesOnNodesPrescribed(self.interface_sub_model_part_origin,
#                                         variable_origin,
#                                         self.vector_values_origin_send)

#         if (self.GiD_output):
#             self.WriteNodalResultsCustom(self.gid_io_origin,
#                                          self.model_part_origin,
#                                          variable_origin,
#                                          output_time)

#         self.mapper.Map(variable_origin,
#                         variable_destination)

#         if (self.GiD_output):
#             self.WriteNodalResultsCustom(self.gid_io_destination,
#                                          self.model_part_destination,
#                                          variable_destination,
#                                          output_time)

#         if (self.set_up_test_2):
#             self.PrintMappedValues(self.interface_sub_model_part_destination,
#                                    variable_destination,
#                                    "Destination_receive Vector")
#         else:
#             self.CheckValuesPrescribed(self.interface_sub_model_part_destination,
#                                        variable_destination,
#                                        self.vector_values_destination_receive)

#     def TestInverseMapNonConstantVectorValues(self, output_time):
#         variable_origin = VELOCITY
#         variable_destination = FORCE

#         self.SetValuesOnNodesPrescribed(self.interface_sub_model_part_destination,
#                                         variable_destination,
#                                         self.vector_values_destination_send)

#         if (self.GiD_output):
#             self.WriteNodalResultsCustom(self.gid_io_destination,
#                                          self.model_part_destination,
#                                          variable_destination,
#                                          output_time)

#         self.mapper.InverseMap(variable_origin,
#                                variable_destination)

#         if (self.GiD_output):
#             self.WriteNodalResultsCustom(self.gid_io_origin,
#                                          self.model_part_origin,
#                                          variable_origin,
#                                          output_time)

#         if (self.set_up_test_2):
#             self.PrintMappedValues(self.interface_sub_model_part_origin,
#                                    variable_origin,
#                                    "Origin_receive Vector")
#         else:
#             self.CheckValuesPrescribed(self.interface_sub_model_part_origin,
#                                        variable_origin,
#                                        self.vector_values_origin_receive)


#     ##### Value Checking Functions
#     def ResetValuesModelParts(self):
#         for node in self.model_part_origin.Nodes:
#             for variable in self.variable_list_scalar:
#                 node.SetSolutionStepValue(variable, 0.0)
#             for variable in self.variable_list_vector:
#                 node.SetSolutionStepValue(variable, Vector([0.0, 0.0, 0.0]))

#         for node in self.model_part_destination.Nodes:
#             for variable in self.variable_list_scalar:
#                 node.SetSolutionStepValue(variable, 0.0)
#             for variable in self.variable_list_vector:
#                 node.SetSolutionStepValue(variable, Vector([0.0, 0.0, 0.0]))


#     def SetValuesOnNodes(self, model_part, variable, value):
#         for node in model_part.Nodes:
#             if (self.parallel_execution):
#                 if (node.GetSolutionStepValue(PARTITION_INDEX) == mpi.rank):
#                     self.SetValuesOnNodesExec(node, variable, value)
#             else:
#                 self.SetValuesOnNodesExec(node, variable, value)

#     def SetValuesOnNodesExec(self, node, variable, value):
#         node.SetSolutionStepValue(variable, value)


#     def SetValuesOnNodesPrescribed(self, model_part, variable, nodal_values):
#         for node in model_part.Nodes:
#             if (self.parallel_execution):
#                 if (node.GetSolutionStepValue(PARTITION_INDEX) == mpi.rank):
#                     self.SetValuesOnNodesPrescribedExec(node, variable, nodal_values)
#             else:
#                 self.SetValuesOnNodesPrescribedExec(node, variable, nodal_values)

#     def SetValuesOnNodesPrescribedExec(self, node, variable, nodal_values):
#         nodal_coords = (node.X, node.Y, node.Z)
#         value_to_prescribe = nodal_values[nodal_coords]
#         if isinstance(value_to_prescribe, tuple):
#             value_to_prescribe = Vector(list(value_to_prescribe))
#         node.SetSolutionStepValue(variable, value_to_prescribe)


#     def CheckValues(self, model_part, variable, value_mapped):
#         for node in model_part.Nodes:
#             if (self.parallel_execution):
#                 if (node.GetSolutionStepValue(PARTITION_INDEX) == mpi.rank):
#                     self.CheckValuesExec(node, variable, value_mapped)
#             else:
#                 self.CheckValuesExec(node, variable, value_mapped)

#     def CheckValuesExec(self, node, variable, value_mapped):
#         value_expected = node.GetSolutionStepValue(variable)
#         self.assertAlmostEqualCustom(value_mapped,value_expected)


#     def CheckValuesPrescribed(self, model_part, variable, nodal_values):
#         for node in model_part.Nodes:
#             if (self.parallel_execution):
#                 if (node.GetSolutionStepValue(PARTITION_INDEX) == mpi.rank):
#                     self.CheckValuesPrescribedExec(node, variable, nodal_values)
#             else:
#                 self.CheckValuesPrescribedExec(node, variable, nodal_values)

#     def CheckValuesPrescribedExec(self, node, variable, nodal_values):
#         value_mapped = node.GetSolutionStepValue(variable)
#         nodal_coords = (node.X, node.Y, node.Z)
#         value_expected = nodal_values[nodal_coords]
#         self.assertAlmostEqualCustom(value_mapped,value_expected)


#     def assertAlmostEqualCustom(self, value_mapped, value_expected):
#         if (isinstance(value_mapped, float) or isinstance(value_mapped, int)): # Variable is a scalar
#             self.assertAlmostEqual(value_mapped,value_expected,4)
#         else: # Variable is a vector
#             for i in range(0,3):
#               self.assertAlmostEqual(value_mapped[i],value_expected[i],4)


#     ##### IO related Functions #####
#     def partition_and_read_model_part(self, model_part_name,
#                                       model_part_input_file,
#                                       size_domain, variable_list,
#                                       number_of_partitions):
#         model_part = ModelPart(model_part_name)
#         for variable in variable_list:
#             model_part.AddNodalSolutionStepVariable(variable)

#         if (number_of_partitions > 1):
#             if (mpi.size > 1):
#                 if (mpi.rank == 0):
#                     model_part_io = ReorderConsecutiveModelPartIO(model_part_input_file)

#                     partitioner = MetisDivideHeterogeneousInputProcess(
#                         model_part_io,
#                         number_of_partitions,
#                         size_domain,
#                         0, # verbosity, set to 1 for more detailed output
#                         True)

#                     partitioner.Execute()

#                 mpi.world.barrier()
#                 model_part_input_file = model_part_input_file + "_" + str(mpi.rank)

#         model_part_io = ModelPartIO(model_part_input_file)
#         model_part_io.ReadModelPart(model_part)

#         if (number_of_partitions > 1):
#             MPICommSetup = SetMPICommunicatorProcess(model_part)
#             MPICommSetup.Execute()

#             ParallelFillComm = ParallelFillCommunicator(model_part.GetRootModelPart())
#             ParallelFillComm.Execute()

#         model_part.ProcessInfo.SetValue(DOMAIN_SIZE, size_domain)
#         model_part.SetBufferSize(1)

#         return model_part

#     def InitializeGiD(self, file_name):
#         # Initialize GidIO
#         output_file_origin = "MapperTests_gid_output/output_" + file_name + "_origin"
#         output_file_destination = "MapperTests_gid_output/output_" + file_name + "_destination"

#         if (self.parallel_execution):
#             output_file_origin += "_r" + str(mpi.rank)
#             output_file_destination += "_r" + str(mpi.rank)

#         gid_mode = GiDPostMode.GiD_PostAscii
#         multifile = MultiFileFlag.MultipleFiles
#         deformed_mesh_flag = WriteDeformedMeshFlag.WriteUndeformed
#         write_conditions = WriteConditionsFlag.WriteConditions

#         self.gid_io_origin = GidIO(output_file_origin, gid_mode, multifile,
#                                    deformed_mesh_flag, write_conditions)
#         self.gid_io_destination = GidIO(output_file_destination, gid_mode, multifile,
#                                         deformed_mesh_flag, write_conditions)

#         # Initialize Results Output
#         self.gid_io_origin.InitializeResults(0, self.model_part_origin.GetMesh())
#         self.gid_io_destination.InitializeResults( 0, self.model_part_destination.GetMesh())

#         # Print original meshes
#         self.write_mesh(self.model_part_origin, self.gid_io_origin)
#         self.write_mesh(self.model_part_destination, self.gid_io_destination)

#     def WriteNodalResultsCustom(self, gid_io, model_part, variable, output_time):
#         if (self.parallel_execution):
#             gid_io.WriteNodalResults(variable, model_part.GetCommunicator().LocalMesh().Nodes, output_time, 0)
#         else:
#             gid_io.WriteNodalResults(variable, model_part.Nodes, output_time, 0)

#     ### Function to write meshes ###
#     def write_mesh(self, model_part, gid_io):
#         gid_io.InitializeMesh(0)
#         gid_io.WriteMesh(model_part.GetMesh())
#         gid_io.FinalizeMesh()

#     def FinalizeGiD(self):
#         # Finalize GidIO
#         self.gid_io_origin.FinalizeResults()
#         self.gid_io_destination.FinalizeResults()


#     ##### Test Set Up related Functions
#     def PrintValuesForJson(self):
#         coordinates_origin_x = []
#         coordinates_origin_y = []
#         coordinates_origin_z = []
#         for node in self.interface_sub_model_part_origin.Nodes:
#             coordinates_origin_x.append(node.X)
#             coordinates_origin_y.append(node.Y)
#             coordinates_origin_z.append(node.Z)

#         coordinates_destination_x = []
#         coordinates_destination_y = []
#         coordinates_destination_z = []
#         for node in self.interface_sub_model_part_destination.Nodes:
#             coordinates_destination_x.append(node.X)
#             coordinates_destination_y.append(node.Y)
#             coordinates_destination_z.append(node.Z)

#         print("\n\nCopy into JSON:\n")
#         print("{")
#         print("    \"Origin_Coordinates\": {")
#         print("        \"X\":", coordinates_origin_x, ",")
#         print("        \"Y\":", coordinates_origin_y, ",")
#         print("        \"Z\":", coordinates_origin_z)
#         print("    },")
#         print("    \"Destination_Coordinates\": {")
#         print("        \"X\":", coordinates_destination_x, ",")
#         print("        \"Y\":", coordinates_destination_y, ",")
#         print("        \"Z\":", coordinates_destination_z)
#         print("    },")
#         print("    \"Origin_send\": {")
#         print("        \"Scalar\":", self.CreateMappingValues(1, coordinates_origin_x, coordinates_origin_y, coordinates_origin_z), ",")
#         print("        \"Vector_X\":", self.CreateMappingValues(2, coordinates_origin_x, coordinates_origin_y, coordinates_origin_z), ",")
#         print("        \"Vector_Y\":", self.CreateMappingValues(3, coordinates_origin_x, coordinates_origin_y, coordinates_origin_z), ",")
#         print("        \"Vector_Z\":", self.CreateMappingValues(4, coordinates_origin_x, coordinates_origin_y, coordinates_origin_z))
#         print("    },")
#         print("    \"Destination_send\": {")
#         print("        \"Scalar\":", self.CreateMappingValues(1, coordinates_destination_x, coordinates_destination_y, coordinates_destination_z), ",")
#         print("        \"Vector_X\":", self.CreateMappingValues(2, coordinates_destination_x, coordinates_destination_y, coordinates_destination_z), ",")
#         print("        \"Vector_Y\":", self.CreateMappingValues(3, coordinates_destination_x, coordinates_destination_y, coordinates_destination_z), ",")
#         print("        \"Vector_Z\":", self.CreateMappingValues(4, coordinates_destination_x, coordinates_destination_y, coordinates_destination_z))
#         print("    },")
#         print("    \"Origin_receive\": {")
#         print("        \"Scalar\": [],")
#         print("        \"Vector_X\": [],")
#         print("        \"Vector_Y\": [],")
#         print("        \"Vector_Z\": []")
#         print("    },")
#         print("    \"Destination_receive\": {")
#         print("        \"Scalar\": [],")
#         print("        \"Vector_X\": [],")
#         print("        \"Vector_Y\": [],")
#         print("        \"Vector_Z\": []")
#         print("    }")
#         print("}")

#         err # needed to get the output

#     def CreateMappingValues(self, index, coords_x, coords_y, coords_z):
#         values_precision = 4
#         values = []
#         value = 0.0

#         for i in range(len(coords_x)):
#             values.append(round(self.CalculateMappingValue(index, coords_x[i], coords_y[i], coords_z[i]), values_precision))

#         return values

#     def CalculateMappingValue(self, index, x, y, z):
#         if (index == 1):
#             value = (x + y + z) * 1e2
#         elif (index == 2):
#             value = ((x**2 + y + z)*1e2)**2
#         elif (index == 3):
#             value = (x**2 + y**2 + z**2)*1e2
#         elif (index == 4):
#             value = ((x + y + z)*1e2)**(3)
#         else:
#             raise("ERROR: Wrong index for MappingValues specified")


#         return value

#     def PrintMappedValues(self, model_part, variable, info_string):
#         values_precision = 6 # higher than the precision of the mapped values bcs of "round()"
#         if ("Vector" in info_string): # Vector Variable
#             values_x = []
#             values_y = []
#             values_z = []
#             for node in model_part.Nodes:
#                 value = node.GetSolutionStepValue(variable)
#                 values_x.append(round(value[0], values_precision))
#                 values_y.append(round(value[1], values_precision))
#                 values_z.append(round(value[2], values_precision))
#             print("\n\n Values for:", info_string)
#             print("        \"Vector_X\":", values_x, ",")
#             print("        \"Vector_Y\":", values_y, ",")
#             print("        \"Vector_Z\":", values_z)

#         else: # Scalar Variable
#             values = []
#             for node in model_part.Nodes:
#                 values.append(round(node.GetSolutionStepValue(variable), values_precision))
#             print("\n\n Values for:", info_string)
#             print("        \"Scalar\":", values, ",")

#     def SetPrescribedValues(self):
#         self.scalar_values_origin_send = {}
#         self.scalar_values_destination_send = {}
#         self.scalar_values_origin_receive = {}
#         self.scalar_values_destination_receive = {}

#         self.vector_values_origin_send = {}
#         self.vector_values_destination_send = {}
#         self.vector_values_origin_receive = {}
#         self.vector_values_destination_receive = {}

#         self.AssignPrescribedValues("Origin",
#                                     "send",
#                                     "Scalar",
#                                     self.scalar_values_origin_send)

#         self.AssignPrescribedValues("Origin",
#                                     "send",
#                                     "Vector",
#                                     self.vector_values_origin_send)

#         self.AssignPrescribedValues("Destination",
#                                     "send",
#                                     "Scalar",
#                                     self.scalar_values_destination_send)

#         self.AssignPrescribedValues("Destination",
#                                     "send",
#                                     "Vector",
#                                     self.vector_values_destination_send)

#         self.AssignPrescribedValues("Origin",
#                                     "receive",
#                                     "Scalar",
#                                     self.scalar_values_origin_receive)

#         self.AssignPrescribedValues("Origin",
#                                     "receive",
#                                     "Vector",
#                                     self.vector_values_origin_receive)

#         self.AssignPrescribedValues("Destination",
#                                     "receive",
#                                     "Scalar",
#                                     self.scalar_values_destination_receive)

#         self.AssignPrescribedValues("Destination",
#                                     "receive",
#                                     "Vector",
#                                     self.vector_values_destination_receive)

#     def AssignPrescribedValues(self, side_interface, direction, value_type, dictionary):
#         coordinates_x = self.results[side_interface + "_Coordinates"]["X"].GetVector()
#         coordinates_y = self.results[side_interface + "_Coordinates"]["Y"].GetVector()
#         coordinates_z = self.results[side_interface + "_Coordinates"]["Z"].GetVector()

#         values = 0
#         if (value_type == "Scalar"):
#             values = self.results[side_interface + "_" + direction]["Scalar"]
#         elif(value_type == "Vector"):
#             values = [ self.results[side_interface + "_" + direction]["Vector_X"],
#                        self.results[side_interface + "_" + direction]["Vector_Y"],
#                        self.results[side_interface + "_" + direction]["Vector_Z"] ]
#         else:
#             raise("ERROR: wrong value_type")

#         i = 0
#         value = 0

#         for x in coordinates_x:
#             # Retreive Coordinates
#             coords = (coordinates_x[i],
#                       coordinates_y[i],
#                       coordinates_z[i])

#             # Retreive Values
#             if (value_type == "Scalar"):
#                 try:
#                     value = values[i].GetDouble()
#                 except: # skip the value in case the results are not yet in the file
#                     pass
#             elif(value_type == "Vector"):
#                 try:
#                     value = (values[0][i].GetDouble(),
#                             values[1][i].GetDouble(),
#                             values[2][i].GetDouble())
#                 except: # skip the value in case the results are not yet in the file
#                     pass

#             dictionary.update({coords : value})
#             i += 1
