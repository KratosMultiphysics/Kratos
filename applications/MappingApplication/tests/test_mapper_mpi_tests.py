from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

import KratosMultiphysics
from KratosMultiphysics.mpi import mpi
import KratosMultiphysics.MetisApplication
import KratosMultiphysics.TrilinosApplication
import KratosMultiphysics.MappingApplication as KratosMapping

import KratosMultiphysics.KratosUnittest as KratosUnittest

from base_mapper_tests import BaseMapperTests
from trilinos_import_model_part_utility import TrilinosImportModelPartUtility

class MapperMPITests(BaseMapperTests, KratosUnittest.TestCase):
    @classmethod
    def _ImportModelPart(cls):
        cls.model_part_origin.AddNodalSolutionStepVariable(
            KratosMultiphysics.PARTITION_INDEX)
        cls.model_part_destination.AddNodalSolutionStepVariable(
            KratosMultiphysics.PARTITION_INDEX)

        origin_settings = KratosMultiphysics.Parameters("""{
            "model_import_settings": {
                "input_type": "mdpa",
                "input_filename": \"""" + cls.input_file_origin + """\",
                "partition_in_memory" : true
            },
            "echo_level" : 0
        }""")
        destination_settings = origin_settings.Clone()
        destination_settings["model_import_settings"]["input_filename"].SetString(
            cls.input_file_destination)

        model_part_import_util_origin = TrilinosImportModelPartUtility(
            cls.model_part_origin, origin_settings)
        model_part_import_util_destination = TrilinosImportModelPartUtility(
            cls.model_part_destination, destination_settings)

        model_part_import_util_origin.ImportModelPart()
        model_part_import_util_destination.ImportModelPart()
        model_part_import_util_origin.CreateCommunicators()
        model_part_import_util_destination.CreateCommunicators()

    def _CreateMapper(self, mapper_settings):
        return KratosMapping.MapperFactory.CreateMPIMapper(
            self.model_part_origin,
            self.model_part_destination,
            mapper_settings)


if __name__ == '__main__':
    KratosUnittest.main()
