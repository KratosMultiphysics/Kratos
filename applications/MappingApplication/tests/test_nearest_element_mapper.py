from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

import KratosMultiphysics as KM
import basic_mapper_tests
import blade_mapping_test

class NearestElementBasicTestsLine(basic_mapper_tests.BasicMapperTests):
    @classmethod
    def setUpClass(cls):
        mapper_params = KM.Parameters("""{
            "mapper_type": "nearest_element",
            "interface_submodel_part_origin": "line_tri",
            "interface_submodel_part_destination": "line_quad",
            "echo_level" : 0
        }""")
        super(NearestElementBasicTestsLine, cls).setUpMapper(mapper_params)

class NearestElementBasicTestsLineSwitchedSides(basic_mapper_tests.BasicMapperTests):
    @classmethod
    def setUpClass(cls):
        mapper_params = KM.Parameters("""{
            "mapper_type": "nearest_element",
            "interface_submodel_part_origin": "line_quad",
            "interface_submodel_part_destination": "line_tri",
            "echo_level" : 0
        }""")
        super(NearestElementBasicTestsLineSwitchedSides, cls).setUpMapper(mapper_params, switch_sides=True)

class NearestElementBasicTestsSurface(basic_mapper_tests.BasicMapperTests):
    @classmethod
    def setUpClass(cls):
        mapper_params = KM.Parameters("""{
            "mapper_type": "nearest_element",
            "interface_submodel_part_origin": "surface_tri",
            "interface_submodel_part_destination": "surface_quad",
            "echo_level" : 0
        }""")
        super(NearestElementBasicTestsSurface, cls).setUpMapper(mapper_params)

class NearestElementBasicTestsSurfaceSwitchedSides(basic_mapper_tests.BasicMapperTests):
    @classmethod
    def setUpClass(cls):
        mapper_params = KM.Parameters("""{
            "mapper_type": "nearest_element",
            "interface_submodel_part_origin": "surface_quad",
            "interface_submodel_part_destination": "surface_tri",
            "echo_level" : 0
        }""")
        super(NearestElementBasicTestsSurfaceSwitchedSides, cls).setUpMapper(mapper_params, switch_sides=True)

class NearestElementBasicTestsVolume(basic_mapper_tests.BasicMapperTests):
    @classmethod
    def setUpClass(cls):
        mapper_params = KM.Parameters("""{
            "mapper_type": "nearest_element",
            "interface_submodel_part_origin": "volume_tri",
            "interface_submodel_part_destination": "volume_quad",
            "echo_level" : 0
        }""")
        super(NearestElementBasicTestsVolume, cls).setUpMapper(mapper_params)

class NearestElementBasicTestsVolumeSwitchedSides(basic_mapper_tests.BasicMapperTests):
    @classmethod
    def setUpClass(cls):
        mapper_params = KM.Parameters("""{
            "mapper_type": "nearest_element",
            "interface_submodel_part_origin": "volume_quad",
            "interface_submodel_part_destination": "volume_tri",
            "echo_level" : 0
        }""")
        super(NearestElementBasicTestsVolumeSwitchedSides, cls).setUpMapper(mapper_params, switch_sides=True)

class NearestElementBladeMapping(blade_mapping_test.BladeMappingTests):
    @classmethod
    def setUpClass(cls):
        mapper_params = KM.Parameters("""{
            "mapper_type": "nearest_element",
            "echo_level" : 0
        }""")
        super(NearestElementBladeMapping, cls).setUpMapper(mapper_params)
        cls.print_output = False

if __name__ == '__main__':
    KM.Logger.GetDefaultOutput().SetSeverity(KM.Logger.Severity.WARNING)
    import KratosMultiphysics.KratosUnittest as KratosUnittest
    KratosUnittest.main()
