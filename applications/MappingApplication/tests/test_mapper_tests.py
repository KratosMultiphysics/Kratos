from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

import KratosMultiphysics
import KratosMultiphysics.MappingApplication as KratosMapping

import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.kratos_utilities as kratos_utils

from base_mapper_tests import BaseMapperTests

class MapperTests(BaseMapperTests, KratosUnittest.TestCase):
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
