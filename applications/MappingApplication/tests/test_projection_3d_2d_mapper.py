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
from math import sin, cos
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

        cls.mapper_type = mapper_parameters["mapper_type"].GetString()

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

    def test_Projection3D2DMapper_Map_non_constant_scalar(self):
        SetHistoricalNonUniformSolutionScalar(self.model_part_origin.Nodes, KM.PRESSURE)
        self.mapper.Map(KM.PRESSURE, KM.TEMPERATURE)
        # mapper_test_case.CheckHistoricalNonUniformValues(self.model_part_destination, KM.TEMPERATURE, GetFilePath(self._GetFileName("map_scalar")))

    # def test_Projection3D2DMapper_InverseMap_non_constant_scalar(self):
    #     SetHistoricalNonUniformSolutionScalar(self.model_part_destination.Nodes, KM.TEMPERATURE)
    #     self.mapper.InverseMap(KM.PRESSURE, KM.TEMPERATURE)
    #     mapper_test_case.CheckHistoricalNonUniformValues(self.model_part_origin, KM.PRESSURE, GetFilePath(self._GetFileName("inverse_map_scalar")))

    # def test_Projection3D2DMapper_Map_non_constant_vector(self):
    #     SetHistoricalNonUniformSolutionVector(self.model_part_origin.Nodes, KM.FORCE)
    #     self.mapper.Map(KM.FORCE, KM.VELOCITY)
    #     mapper_test_case.CheckHistoricalNonUniformValues(self.model_part_destination, KM.VELOCITY, GetFilePath(self._GetFileName("map_vector")))

    # def test_Projection3D2DMapper_InverseMap_non_constant_vector(self):
    #     SetHistoricalNonUniformSolutionVector(self.model_part_destination.Nodes, KM.VELOCITY)
    #     self.mapper.InverseMap(KM.FORCE, KM.VELOCITY)
    #     mapper_test_case.CheckHistoricalNonUniformValues(self.model_part_origin, KM.FORCE, GetFilePath(self._GetFileName("inverse_map_vector")))

    def _CheckHistoricalUniformValuesScalar(self, nodes, variable, exp_value):
        for node in nodes:
            self.assertAlmostEqual(node.GetSolutionStepValue(variable), exp_value)

    def _CheckHistoricalUniformValuesVector(self, nodes, variable, exp_value):
        for node in nodes:
            self.assertVectorAlmostEqual(node.GetSolutionStepValue(variable), exp_value)

    def _CheckUniformValuesScalar(self, entities, variable, exp_value):
        for entity in entities:
            self.assertAlmostEqual(entity.GetValue(variable), exp_value)

    def _CheckUniformValuesVector(self, entities, variable, exp_value):
        for entity in entities:
            self.assertVectorAlmostEqual(entity.GetValue(variable), exp_value)

    def _GetFileName(self, file_appendix):
        return os.path.join("result_files", self.mapper_type, self.__class__.__name__ + "_" + file_appendix)

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

def GetNodes(model_part):
    return model_part.Nodes

if __name__ == '__main__':
    KM.Logger.GetDefaultOutput().SetSeverity(KM.Logger.Severity.WARNING)
    import KratosMultiphysics.KratosUnittest as KratosUnittest
    KratosUnittest.main()
