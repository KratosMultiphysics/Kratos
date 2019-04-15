from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

import KratosMultiphysics as KM
import KratosMultiphysics.MappingApplication as KratosMapping
import basic_mapper_tests_new

class NearestElementBasicTestsLine(basic_mapper_tests_new.BasicMapperTests):
    @classmethod
    def setUpClass(cls):
        mapper_params = KM.Parameters("""{
            "mapper_type": "nearest_element",
            "interface_submodel_part_origin": "LineLoad3D_mapping_line_tri",
            "interface_submodel_part_destination": "LineLoad3D_mapping_line_quad",
            "echo_level" : 0
        }""")
        super(NearestElementBasicTestsLine, cls).setUpClass(mapper_params)

class NearestElementBasicTestsLineSwitchedSides(basic_mapper_tests_new.BasicMapperTests):
    @classmethod
    def setUpClass(cls):
        mapper_params = KM.Parameters("""{
            "mapper_type": "nearest_element",
            "interface_submodel_part_origin": "LineLoad3D_mapping_line_quad",
            "interface_submodel_part_destination": "LineLoad3D_mapping_line_tri",
            "echo_level" : 0
        }""")
        super(NearestElementBasicTestsLineSwitchedSides, cls).setUpClass(mapper_params, switch_sides=True)

class NearestElementBasicTestsSurface(basic_mapper_tests_new.BasicMapperTests):
    @classmethod
    def setUpClass(cls):
        mapper_params = KM.Parameters("""{
            "mapper_type": "nearest_element",
            "interface_submodel_part_origin": "SurfaceLoad3D_mapping_surface_tri",
            "interface_submodel_part_destination": "SurfaceLoad3D_mapping_surface_quad",
            "echo_level" : 0
        }""")
        super(NearestElementBasicTestsSurface, cls).setUpClass(mapper_params)

class NearestElementBasicTestsSurfaceSwitchedSides(basic_mapper_tests_new.BasicMapperTests):
    @classmethod
    def setUpClass(cls):
        mapper_params = KM.Parameters("""{
            "mapper_type": "nearest_element",
            "interface_submodel_part_origin": "SurfaceLoad3D_mapping_surface_quad",
            "interface_submodel_part_destination": "SurfaceLoad3D_mapping_surface_tri",
            "echo_level" : 0
        }""")
        super(NearestElementBasicTestsSurfaceSwitchedSides, cls).setUpClass(mapper_params, switch_sides=True)

class NearestElementBasicTestsVolume(basic_mapper_tests_new.BasicMapperTests):
    @classmethod
    def setUpClass(cls):
        mapper_params = KM.Parameters("""{
            "mapper_type": "nearest_element",
            "interface_submodel_part_origin": "Parts_domain_tri",
            "interface_submodel_part_destination": "Parts_domain_quad",
            "echo_level" : 0
        }""")
        super(NearestElementBasicTestsVolume, cls).setUpClass(mapper_params)

class NearestElementBasicTestsVolumeSwitchedSides(basic_mapper_tests_new.BasicMapperTests):
    @classmethod
    def setUpClass(cls):
        mapper_params = KM.Parameters("""{
            "mapper_type": "nearest_element",
            "interface_submodel_part_origin": "Parts_domain_quad",
            "interface_submodel_part_destination": "Parts_domain_tri",
            "echo_level" : 0
        }""")
        super(NearestElementBasicTestsVolumeSwitchedSides, cls).setUpClass(mapper_params, switch_sides=True)

if __name__ == '__main__':
    import sys
    if "--using-mpi" in sys.argv:
        # this is a hack until unittest supports MPI
        from KratosMultiphysics import mpi
        sys.argv.remove("--using-mpi")
    KM.Logger.GetDefaultOutput().SetSeverity(KM.Logger.Severity.WARNING)
    import KratosMultiphysics.KratosUnittest as KratosUnittest
    KratosUnittest.main()
