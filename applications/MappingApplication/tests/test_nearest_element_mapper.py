import KratosMultiphysics as KM
import KratosMultiphysics.MappingApplication # registering the mappers
import basic_mapper_tests
import blade_mapping_test

class BasicTestsLine(basic_mapper_tests.BasicMapperTests):
    @classmethod
    def setUpClass(cls):
        mapper_params = KM.Parameters("""{
            "mapper_type": "nearest_element",
            "interface_submodel_part_origin": "line_tri",
            "interface_submodel_part_destination": "line_quad",
            "echo_level" : 0
        }""")
        super().setUpMapper(mapper_params)

class BasicTestsLineInitialConfig(basic_mapper_tests.BasicMapperTests):
    @classmethod
    def setUpClass(cls):
        mapper_params = KM.Parameters("""{
            "mapper_type": "nearest_element",
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
            "mapper_type": "nearest_element",
            "interface_submodel_part_origin": "line_quad",
            "interface_submodel_part_destination": "line_tri",
            "echo_level" : 0
        }""")
        super().setUpMapper(mapper_params, switch_sides=True)

class BasicTestsSurface(basic_mapper_tests.BasicMapperTests):
    @classmethod
    def setUpClass(cls):
        mapper_params = KM.Parameters("""{
            "mapper_type": "nearest_element",
            "interface_submodel_part_origin": "surface_tri",
            "interface_submodel_part_destination": "surface_quad",
            "echo_level" : 0
        }""")
        super().setUpMapper(mapper_params)

class BasicTestsSurfaceSwitchedSides(basic_mapper_tests.BasicMapperTests):
    @classmethod
    def setUpClass(cls):
        mapper_params = KM.Parameters("""{
            "mapper_type": "nearest_element",
            "interface_submodel_part_origin": "surface_quad",
            "interface_submodel_part_destination": "surface_tri",
            "echo_level" : 0
        }""")
        super().setUpMapper(mapper_params, switch_sides=True)

class BasicTestsVolume(basic_mapper_tests.BasicMapperTests):
    @classmethod
    def setUpClass(cls):
        mapper_params = KM.Parameters("""{
            "mapper_type": "nearest_element",
            "interface_submodel_part_origin": "volume_tri",
            "interface_submodel_part_destination": "volume_quad",
            "echo_level" : 0
        }""")
        super().setUpMapper(mapper_params)

class BasicTestsVolumeSwitchedSides(basic_mapper_tests.BasicMapperTests):
    @classmethod
    def setUpClass(cls):
        mapper_params = KM.Parameters("""{
            "mapper_type": "nearest_element",
            "interface_submodel_part_origin": "volume_quad",
            "interface_submodel_part_destination": "volume_tri",
            "echo_level" : 0
        }""")
        super().setUpMapper(mapper_params, switch_sides=True)

class BladeMapping(blade_mapping_test.BladeMappingTests):
    @classmethod
    def setUpClass(cls):
        mapper_params = KM.Parameters("""{
            "mapper_type": "nearest_element",
            "echo_level" : 0
        }""")
        super().setUpMapper(mapper_params)
        cls.print_output = False

class BladeMappingSerialModelPart(blade_mapping_test.BladeMappingTestsSerialModelPart):
    @classmethod
    def setUpClass(cls):
        mapper_params = KM.Parameters("""{
            "mapper_type": "nearest_element",
            "echo_level" : 0
        }""")
        super().setUpMapper(mapper_params)
        cls.print_output = False

class BladeMappingAllRanksExceptLast(blade_mapping_test.BladeMappingTestsAllRanksExceptLast):
    @classmethod
    def setUpClass(cls):
        mapper_params = KM.Parameters("""{
            "mapper_type": "nearest_element",
            "echo_level" : 0
        }""")
        super().setUpMapper(mapper_params)
        cls.print_output = False

class BladeMappingAllRanksExceptFirst(blade_mapping_test.BladeMappingTestsAllRanksExceptFirst):
    @classmethod
    def setUpClass(cls):
        mapper_params = KM.Parameters("""{
            "mapper_type": "nearest_element",
            "echo_level" : 0
        }""")
        super().setUpMapper(mapper_params)
        cls.print_output = False

class BladeMappingUnevenRanks(blade_mapping_test.BladeMappingTestsUnevenRanks):
    @classmethod
    def setUpClass(cls):
        mapper_params = KM.Parameters("""{
            "mapper_type": "nearest_element",
            "echo_level" : 0
        }""")
        super().setUpMapper(mapper_params)
        cls.print_output = False

if __name__ == '__main__':
    KM.Logger.GetDefaultOutput().SetSeverity(KM.Logger.Severity.WARNING)
    import KratosMultiphysics.KratosUnittest as KratosUnittest
    KratosUnittest.main()
