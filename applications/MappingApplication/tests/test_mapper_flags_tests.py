from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

import KratosMultiphysics
import KratosMultiphysics.MappingApplication as KratosMapping

import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.kratos_utilities as kratos_utils

from base_mapper_tests import MapperTestsBase

class MapperFlagsTests(MapperTestsBase, KratosUnittest.TestCase):

    def test_SWAP_SIGN(self):
        mapper_settings = KratosMultiphysics.Parameters("""{
            "mapper_type": "nearest_neighbor",
            "interface_submodel_part_origin": "LineLoad3D_mapping_line_tri",
            "interface_submodel_part_destination": "LineLoad3D_mapping_line_quad"
        }""")
        values_file_name = "nearest_neighbor_line"

    def test_ADD_VALUES(self):
        mapper_settings = KratosMultiphysics.Parameters("""{
            "mapper_type": "nearest_neighbor",
            "interface_submodel_part_origin": "LineLoad3D_mapping_line_tri",
            "interface_submodel_part_destination": "LineLoad3D_mapping_line_quad"
        }""")
        values_file_name = "nearest_neighbor_line"

    def test_TO_NON_HISTORICAL(self):
        mapper_settings = KratosMultiphysics.Parameters("""{
            "mapper_type": "nearest_neighbor",
            "interface_submodel_part_origin": "LineLoad3D_mapping_line_tri",
            "interface_submodel_part_destination": "LineLoad3D_mapping_line_quad"
        }""")
        values_file_name = "nearest_neighbor_line"

        # self.__ExecuteMapperTests(mapper_settings, values_file_name)

    def test_FROM_NON_HISTORICAL(self):
        mapper_settings = KratosMultiphysics.Parameters("""{
            "mapper_type": "nearest_neighbor",
            "interface_submodel_part_origin": "LineLoad3D_mapping_line_tri",
            "interface_submodel_part_destination": "LineLoad3D_mapping_line_quad"
        }""")
        values_file_name = "nearest_neighbor_line"

        # self.__ExecuteMapperTests(mapper_settings, values_file_name)

    def test_NON_HISTORICAL(self):
        mapper_settings = KratosMultiphysics.Parameters("""{
            "mapper_type": "nearest_neighbor",
            "interface_submodel_part_origin": "LineLoad3D_mapping_line_tri",
            "interface_submodel_part_destination": "LineLoad3D_mapping_line_quad"
        }""")
        values_file_name = "nearest_neighbor_line"

        # self.__ExecuteMapperTests(mapper_settings, values_file_name)

    def test_SWAP_SIGN_and_ADD_VALUES(self):
        mapper_settings = KratosMultiphysics.Parameters("""{
            "mapper_type": "nearest_neighbor",
            "interface_submodel_part_origin": "LineLoad3D_mapping_line_tri",
            "interface_submodel_part_destination": "LineLoad3D_mapping_line_quad"
        }""")
        values_file_name = "nearest_neighbor_line"

        # self.__ExecuteMapperTests(mapper_settings, values_file_name)

    def test_SWAP_SIGN_and_ADD_VALUES_and_CONSERVATIVE(self):
        mapper_settings = KratosMultiphysics.Parameters("""{
            "mapper_type": "nearest_neighbor",
            "interface_submodel_part_origin": "LineLoad3D_mapping_line_tri",
            "interface_submodel_part_destination": "LineLoad3D_mapping_line_quad"
        }""")
        values_file_name = "nearest_neighbor_line"

        # self.__ExecuteMapperTests(mapper_settings, values_file_name)

    def _CreateMapper(self, mapper_settings):
        '''In an MPI-test this function returns a distributed mapper
        In the base-class it returns a regular Mapper
        '''
        return KratosMapping.MapperFactory.CreateMapper(
            self.model_part_origin,
            self.model_part_destination,
            mapper_settings)


if __name__ == '__main__':
    KratosUnittest.main()
