import KratosMultiphysics as KM
data_comm = KM.Testing.GetDefaultDataCommunicator()
if data_comm.IsDistributed():
    from KratosMultiphysics.MappingApplication import MPIExtension as MappingMPIExtension

import mapper_test_case
from math import sin
import os

def GetFilePath(file_name):
    return os.path.join(os.path.dirname(os.path.realpath(__file__)), file_name)

class QuadraticMapperTests(mapper_test_case.MapperTestCase):
    '''This class contains the basic tests for the quadratic mappers.
    '''

    @classmethod
    def setUpMapper(cls, mapper_parameters, switch_sides=False):
        if switch_sides:
            super().setUpModelParts("cube_quadratic", "cube_linear")
        else:
            super().setUpModelParts("cube_linear", "cube_quadratic")
        # TODO ATTENTION: currently the MapperFactory removes some keys, hence those checks have to be done beforehand => improve this!

        cls.mapper_type = mapper_parameters["mapper_type"].GetString()
        cls.mapper_parameters = mapper_parameters.Clone()

        if mapper_parameters.Has("interface_submodel_part_origin"):
            cls.interface_model_part_origin = cls.model_part_origin.GetSubModelPart(
                mapper_parameters["interface_submodel_part_origin"].GetString())
        else:
            cls.interface_model_part_origin = cls.model_part_origin

        if mapper_parameters.Has("interface_submodel_part_destination"):
            cls.interface_model_part_destination = cls.model_part_destination.GetSubModelPart(
                mapper_parameters["interface_submodel_part_destination"].GetString())
        else:
            cls.interface_model_part_destination = cls.model_part_destination

        if data_comm.IsDistributed():
            cls.mapper = MappingMPIExtension.MPIMapperFactory.CreateMapper(
                cls.model_part_origin, cls.model_part_destination, mapper_parameters)
        else:
            cls.mapper = KM.MapperFactory.CreateMapper(
                cls.model_part_origin, cls.model_part_destination, mapper_parameters)

    def test_quadratic_Map_non_constant_scalar(self):
        SetHistoricalNonUniformSolutionScalar(self.interface_model_part_origin.Nodes, KM.PRESSURE)
        self.mapper.Map(KM.PRESSURE, KM.TEMPERATURE)
        mapper_test_case.CheckHistoricalNonUniformValues(self.interface_model_part_destination, KM.TEMPERATURE, GetFilePath(self._GetFileName("map_scalar_quadratic")), output_reference_solution=False)

    def _GetFileName(self, file_appendix):
        return os.path.join("result_files", self.mapper_type, self.__class__.__name__ + "_" + file_appendix)

def SetHistoricalNonUniformSolutionScalar(nodes, variable):
    for node in nodes:
        val = 12*sin(node.X0) + node.Y0*15 + 22*node.Z0
        node.SetSolutionStepValue(variable, val)

def GetNodes(model_part):
    return model_part.Nodes
