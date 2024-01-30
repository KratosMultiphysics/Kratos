import KratosMultiphysics as KM
import KratosMultiphysics.MappingApplication # registering the mappers
import KratosMultiphysics.KratosUnittest as KratosUnittest
from KratosMultiphysics import ParallelEnvironment, IsDistributedRun
import basic_mapper_tests
import blade_mapping_test
data_comm = KM.Testing.GetDefaultDataCommunicator()
if data_comm.IsDistributed():
    from KratosMultiphysics.MappingApplication import MPIExtension as MappingMPIExtension

# Additional imports for corner cases
import mapper_test_case
import os
from math import sin, cos

def GetFilePath(file_name):
    return os.path.join(os.path.dirname(os.path.realpath(__file__)), file_name)

class BasicTestsLine(basic_mapper_tests.BasicMapperTests):
    @classmethod
    def setUpClass(cls):
        mapper_params = KM.Parameters("""{
            "mapper_type": "nearest_neighbor",
            "interface_submodel_part_origin": "line_tri",
            "interface_submodel_part_destination": "line_quad",
            "echo_level" : 0
        }""")
        super().setUpMapper(mapper_params)

class BasicTestsLineInitialConfig(basic_mapper_tests.BasicMapperTests):
    @classmethod
    def setUpClass(cls):
        mapper_params = KM.Parameters("""{
            "mapper_type": "nearest_neighbor",
            "interface_submodel_part_origin": "line_tri",
            "interface_submodel_part_destination": "line_quad",
            "use_initial_configuration" : true,
            "echo_level" : 0
        }""")
        super().setUpMapper(mapper_params)
        for node in cls.model_part_origin.Nodes:
            node.X = node.X + 111.1
            node.Y = node.Y - 693.1
            node.Z = node.Z + 15698

    def _GetFileName(self, file_appendix):
        file_name = super()._GetFileName(file_appendix)
        return file_name.replace("InitialConfig", "")

class BasicTestsLineSwitchedSides(basic_mapper_tests.BasicMapperTests):
    @classmethod
    def setUpClass(cls):
        mapper_params = KM.Parameters("""{
            "mapper_type": "nearest_neighbor",
            "interface_submodel_part_origin": "line_quad",
            "interface_submodel_part_destination": "line_tri",
            "echo_level" : 0
        }""")
        super().setUpMapper(mapper_params, switch_sides=True)

class BasicTestsSurface(basic_mapper_tests.BasicMapperTests):
    @classmethod
    def setUpClass(cls):
        mapper_params = KM.Parameters("""{
            "mapper_type": "nearest_neighbor",
            "interface_submodel_part_origin": "surface_tri",
            "interface_submodel_part_destination": "surface_quad",
            "echo_level" : 0
        }""")
        super().setUpMapper(mapper_params)

class BasicTestsSurfaceSwitchedSides(basic_mapper_tests.BasicMapperTests):
    @classmethod
    def setUpClass(cls):
        mapper_params = KM.Parameters("""{
            "mapper_type": "nearest_neighbor",
            "interface_submodel_part_origin": "surface_quad",
            "interface_submodel_part_destination": "surface_tri",
            "echo_level" : 0
        }""")
        super().setUpMapper(mapper_params, switch_sides=True)

class BasicTestsVolume(basic_mapper_tests.BasicMapperTests):
    @classmethod
    def setUpClass(cls):
        mapper_params = KM.Parameters("""{
            "mapper_type": "nearest_neighbor",
            "echo_level" : 0
        }""")
        super().setUpMapper(mapper_params)

    def test_Is_not_conforming(self):
        non_conform_parameters = KM.Parameters("""{
            "mapper_type": "nearest_neighbor",
            "echo_level" : 0,
            "search_settings" : {
                "search_radius": 1e-8,
                "max_num_search_iterations": 2
            }
        }""")

        if data_comm.IsDistributed():
            map_creator = MappingMPIExtension.MPIMapperFactory.CreateMapper
        else:
            map_creator = KM.MapperFactory.CreateMapper

        non_conform_mapper = map_creator(
            self.model_part_origin,
            self.model_part_destination,
            non_conform_parameters
        )

        is_conforming = non_conform_mapper.AreMeshesConforming()
        self.assertFalse(is_conforming)

    def test_Is_conforming(self):
        conform_parameters = KM.Parameters("""{
            "mapper_type": "nearest_neighbor",
            "echo_level" : 0
        }""")

        if data_comm.IsDistributed():
            map_creator = MappingMPIExtension.MPIMapperFactory.CreateMapper
        else:
            map_creator = KM.MapperFactory.CreateMapper

        non_conform_mapper = map_creator(
            self.model_part_origin,
            self.model_part_origin,
            conform_parameters
        )

        is_conforming = non_conform_mapper.AreMeshesConforming()
        self.assertTrue(is_conforming)

class BasicTestsVolumeSwitchedSides(basic_mapper_tests.BasicMapperTests):
    @classmethod
    def setUpClass(cls):
        mapper_params = KM.Parameters("""{
            "mapper_type": "nearest_neighbor",
            "echo_level" : 0
        }""")
        super().setUpMapper(mapper_params, switch_sides=True)

class BladeMapping(blade_mapping_test.BladeMappingTests):
    @classmethod
    def setUpClass(cls):
        mapper_params = KM.Parameters("""{
            "mapper_type": "nearest_neighbor",
            "echo_level" : 0
        }""")
        super().setUpMapper(mapper_params)
        cls.print_output = False

class BladeMappingSerialModelPart(blade_mapping_test.BladeMappingTestsSerialModelPart):
    @classmethod
    def setUpClass(cls):
        mapper_params = KM.Parameters("""{
            "mapper_type": "nearest_neighbor",
            "echo_level" : 0
        }""")
        super().setUpMapper(mapper_params)
        cls.print_output = False

class BladeMappingAllRanksExceptLast(blade_mapping_test.BladeMappingTestsAllRanksExceptLast):
    @classmethod
    def setUpClass(cls):
        mapper_params = KM.Parameters("""{
            "mapper_type": "nearest_neighbor",
            "echo_level" : 0
        }""")
        super().setUpMapper(mapper_params)
        cls.print_output = False

class BladeMappingAllRanksExceptFirst(blade_mapping_test.BladeMappingTestsAllRanksExceptFirst):
    @classmethod
    def setUpClass(cls):
        mapper_params = KM.Parameters("""{
            "mapper_type": "nearest_neighbor",
            "echo_level" : 0
        }""")
        super().setUpMapper(mapper_params)
        cls.print_output = False

class BladeMappingUnevenRanks(blade_mapping_test.BladeMappingTestsUnevenRanks):
    @classmethod
    def setUpClass(cls):
        mapper_params = KM.Parameters("""{
            "mapper_type": "nearest_neighbor",
            "echo_level" : 0
        }""")
        super().setUpMapper(mapper_params)
        cls.print_output = False

class CornerCaseNearestNeighbor(mapper_test_case.MapperTestCase):
    '''This class contains the tests for corner case for the NearestNeighbor mapper
    '''

    @classmethod
    def setUpMapper(cls, mapper_parameters):
        mdpa_1 = "3D_blocks_mesh1"
        mdpa_2 = "3D_blocks_mesh2"
        super().setUpModelParts(mdpa_1, mdpa_2)

        # Will be 2D
        cls.model_part_origin = cls.model_part_origin.GetSubModelPart("Parts_2D")
        cls.model_part_destination = cls.model_part_destination.GetSubModelPart("Parts_2D")

        cls.mapper_type = mapper_parameters["mapper_type"].GetString()

        if IsDistributedRun():
            cls.mapper = MappingMPIExtension.MPIMapperFactory.CreateMapper(cls.model_part_origin, cls.model_part_destination, mapper_parameters)
        else:
            cls.mapper = KM.MapperFactory.CreateMapper(cls.model_part_origin, cls.model_part_destination, mapper_parameters)

    @classmethod
    def setUpClass(cls):
        mapper_params = KM.Parameters("""{
            "mapper_type"                        : "nearest_neighbor",
            "echo_level"                         : 0
        }""")
        cls.setUpMapper(mapper_params)

    def _GetFileName(self, file_appendix):
        return os.path.join("result_files", self.mapper_type, self.__class__.__name__ + "_" + file_appendix)

    def test_CornerCaseNearestNeighbor_Map_non_constant_scalar(self):
        SetHistoricalNonUniformSolutionScalar(self.model_part_origin.Nodes, KM.PRESSURE)
        self.mapper.Map(KM.PRESSURE, KM.TEMPERATURE)
        mapper_test_case.CheckHistoricalNonUniformValues(self.model_part_destination, KM.TEMPERATURE, GetFilePath(self._GetFileName("map_scalar")))
        #mapper_test_case.VtkOutputNodesHistorical(self.model_part_destination, KM.TEMPERATURE)

    def test_CornerCaseNearestNeighbor_InverseMap_non_constant_scalar(self):
        SetHistoricalNonUniformSolutionScalar(self.model_part_destination.Nodes, KM.TEMPERATURE)
        self.mapper.InverseMap(KM.PRESSURE, KM.TEMPERATURE)
        mapper_test_case.CheckHistoricalNonUniformValues(self.model_part_origin, KM.PRESSURE, GetFilePath(self._GetFileName("inverse_map_scalar")))
        #mapper_test_case.VtkOutputNodesHistorical(self.model_part_origin, KM.PRESSURE)

def SetHistoricalNonUniformSolutionScalar(nodes, variable):
    for node in nodes:
        val = 12*sin(node.X0) + node.Y0*15
        node.SetSolutionStepValue(variable, val)

def SetHistoricalNonUniformSolutionVector(nodes, variable):
    for node in nodes:
        val_1 = 12*sin(node.X0) + node.Y0*15
        val_2 = 33*cos(node.X0) + node.Y0*5
        val_3 = 12*sin(node.Y0) + 22*node.X0
        node.SetSolutionStepValue(variable, KM.Vector([val_1, val_2, val_3]))

if __name__ == '__main__':
    KM.Logger.GetDefaultOutput().SetSeverity(KM.Logger.Severity.WARNING)
    import KratosMultiphysics.KratosUnittest as KratosUnittest
    KratosUnittest.main()
