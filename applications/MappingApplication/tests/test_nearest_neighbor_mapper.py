from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

import KratosMultiphysics as KM
import basic_mapper_tests
import blade_mapping_test

class NearestNeighborBasicTestsLine(basic_mapper_tests.BasicMapperTests):
    @classmethod
    def setUpClass(cls):
        mapper_params = KM.Parameters("""{
            "mapper_type": "nearest_neighbor",
            "interface_submodel_part_origin": "line_tri",
            "interface_submodel_part_destination": "line_quad",
            "echo_level" : 0
        }""")
        super(NearestNeighborBasicTestsLine, cls).setUpMapper(mapper_params)

class NearestNeighborBasicTestsLineSwitchedSides(basic_mapper_tests.BasicMapperTests):
    @classmethod
    def setUpClass(cls):
        mapper_params = KM.Parameters("""{
            "mapper_type": "nearest_neighbor",
            "interface_submodel_part_origin": "line_quad",
            "interface_submodel_part_destination": "line_tri",
            "echo_level" : 0
        }""")
        super(NearestNeighborBasicTestsLineSwitchedSides, cls).setUpMapper(mapper_params, switch_sides=True)

class NearestNeighborBasicTestsSurface(basic_mapper_tests.BasicMapperTests):
    @classmethod
    def setUpClass(cls):
        mapper_params = KM.Parameters("""{
            "mapper_type": "nearest_neighbor",
            "interface_submodel_part_origin": "surface_tri",
            "interface_submodel_part_destination": "surface_quad",
            "echo_level" : 0
        }""")
        super(NearestNeighborBasicTestsSurface, cls).setUpMapper(mapper_params)

class NearestNeighborBasicTestsSurfaceSwitchedSides(basic_mapper_tests.BasicMapperTests):
    @classmethod
    def setUpClass(cls):
        mapper_params = KM.Parameters("""{
            "mapper_type": "nearest_neighbor",
            "interface_submodel_part_origin": "surface_quad",
            "interface_submodel_part_destination": "surface_tri",
            "echo_level" : 0
        }""")
        super(NearestNeighborBasicTestsSurfaceSwitchedSides, cls).setUpMapper(mapper_params, switch_sides=True)

class NearestNeighborBasicTestsVolume(basic_mapper_tests.BasicMapperTests):
    @classmethod
    def setUpClass(cls):
        mapper_params = KM.Parameters("""{
            "mapper_type": "nearest_neighbor",
            "echo_level" : 0
        }""")
        super(NearestNeighborBasicTestsVolume, cls).setUpMapper(mapper_params)

class NearestNeighborBasicTestsVolumeSwitchedSides(basic_mapper_tests.BasicMapperTests):
    @classmethod
    def setUpClass(cls):
        mapper_params = KM.Parameters("""{
            "mapper_type": "nearest_neighbor",
            "echo_level" : 0
        }""")
        super(NearestNeighborBasicTestsVolumeSwitchedSides, cls).setUpMapper(mapper_params, switch_sides=True)

class NearestNeighborBladeMapping(blade_mapping_test.BladeMappingTests):
    @classmethod
    def setUpClass(cls):
        mapper_params = KM.Parameters("""{
            "mapper_type": "nearest_neighbor",
            "echo_level" : 0
        }""")
        super(NearestNeighborBladeMapping, cls).setUpMapper(mapper_params)
        cls.print_output = False

if __name__ == '__main__':
    KM.Logger.GetDefaultOutput().SetSeverity(KM.Logger.Severity.WARNING)
    import KratosMultiphysics.KratosUnittest as KratosUnittest
    KratosUnittest.main()
