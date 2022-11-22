import KratosMultiphysics as KM
import KratosMultiphysics.MappingApplication # registering the mappers
from KratosMultiphysics import KratosUnittest
default_data_comm = KM.Testing.GetDefaultDataCommunicator()
if default_data_comm.IsDistributed():
    from KratosMultiphysics import mpi as KratosMPI
    from KratosMultiphysics.MappingApplication import MPIExtension as MappingMPIExtension

from KratosMultiphysics.testing import utilities as testing_utils
import mapper_test_case
import os
from sys import version_info as py_version_info

def GetFilePath(file_name):
    return os.path.join(os.path.dirname(os.path.realpath(__file__)), file_name)

class Projection3D2DMapper(mapper_test_case.MapperTestCase):
    '''This class contains the tests for 3D-2D projections
    '''

    @classmethod
    def setUpMapper(cls, mapper_parameters):
        mdpa_1 = "3D_blocks_mesh1"
        mdpa_2 = "3D_blocks_mesh2"
        super().setUpModelParts(mdpa_1, mdpa_2)

        if default_data_comm.IsDistributed():
            map_creator = MappingMPIExtension.MPIMapperFactory.CreateMapper
        else:
            map_creator = KM.MapperFactory.CreateMapper

        cls.mapper = map_creator(
            cls.model_part_origin,
            cls.model_part_destination,
            mapper_parameters.Clone()
            )

    @classmethod
    def setUpClass(cls):
        mapper_params = KM.Parameters("""{
            "mapper_type" : "projection_3D_2D",
            "base_mapper" : "nearest_neighbor",
            "echo_level"  : 0
        }""")
        cls.setUpMapper(mapper_params)

    def test_Projection3D2DMapper(self):
        pass

if __name__ == '__main__':
    KM.Logger.GetDefaultOutput().SetSeverity(KM.Logger.Severity.WARNING)
    import KratosMultiphysics.KratosUnittest as KratosUnittest
    KratosUnittest.main()
