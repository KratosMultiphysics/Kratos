import KratosMultiphysics as KM
from KratosMultiphysics import ParallelEnvironment, IsDistributedRun
import KratosMultiphysics.MappingApplication # registering the mappers
import KratosMultiphysics.KratosUnittest as KratosUnittest
if IsDistributedRun():
    from KratosMultiphysics.MappingApplication import MPIExtension as MappingMPIExtension

import mapper_test_case
import os
from math import sin, cos

def GetFilePath(file_name):
    return os.path.join(os.path.dirname(os.path.realpath(__file__)), file_name)

class Projection3D2DMapperNearestNeighborOrigin2D(mapper_test_case.MapperTestCase):
    '''This class contains the tests for 3D-2D projections for the NearestNeighbor mapper (origin 2D)
    '''

    @classmethod
    def setUpMapper(cls, mapper_parameters):
        mdpa_1 = "3D_blocks_mesh1"
        mdpa_2 = "3D_blocks_mesh2"
        super().setUpModelParts(mdpa_1, mdpa_2)

        # Origin will be 2D
        cls.model_part_origin = cls.model_part_origin.GetSubModelPart("Parts_2D")

        cls.mapper_type = mapper_parameters["mapper_type"].GetString()

        if IsDistributedRun():
            cls.mapper = MappingMPIExtension.MPIMapperFactory.CreateMapper(cls.model_part_origin, cls.model_part_destination, mapper_parameters)
        else:
            cls.mapper = KM.MapperFactory.CreateMapper(cls.model_part_origin, cls.model_part_destination, mapper_parameters)

    @classmethod
    def setUpClass(cls):
        mapper_params = KM.Parameters("""{
            "mapper_type"                        : "projection_3D_2D",
            "base_mapper"                        : "nearest_neighbor",
            "echo_level"                         : 0
        }""")
        cls.setUpMapper(mapper_params)

    def _GetFileName(self, file_appendix):
        return os.path.join("result_files", self.mapper_type, self.__class__.__name__ + "_" + file_appendix)

    def test_Projection3D2DMapper_Map_non_constant_scalar(self):
        SetHistoricalNonUniformSolutionScalar(self.model_part_origin.Nodes, KM.PRESSURE)
        self.mapper.Map(KM.PRESSURE, KM.TEMPERATURE)
        mapper_test_case.CheckHistoricalNonUniformValues(self.model_part_destination, KM.TEMPERATURE, GetFilePath(self._GetFileName("map_scalar")))
        #mapper_test_case.VtkOutputNodesHistorical(self.model_part_destination, KM.TEMPERATURE)

    def test_Projection3D2DMapper_InverseMap_non_constant_scalar(self):
        SetHistoricalNonUniformSolutionScalar(self.model_part_destination.Nodes, KM.TEMPERATURE)
        self.mapper.InverseMap(KM.PRESSURE, KM.TEMPERATURE)
        mapper_test_case.CheckHistoricalNonUniformValues(self.model_part_origin, KM.PRESSURE, GetFilePath(self._GetFileName("inverse_map_scalar")))
        #mapper_test_case.VtkOutputNodesHistorical(self.model_part_origin, KM.PRESSURE)

    def test_Projection3D2DMapper_Map_non_constant_vector(self):
        SetHistoricalNonUniformSolutionVector(self.model_part_origin.Nodes, KM.FORCE)
        self.mapper.Map(KM.FORCE, KM.VELOCITY)
        mapper_test_case.CheckHistoricalNonUniformValues(self.model_part_destination, KM.VELOCITY, GetFilePath(self._GetFileName("map_vector")))
        #mapper_test_case.VtkOutputNodesHistorical(self.model_part_destination, KM.VELOCITY)

    def test_Projection3D2DMapper_InverseMap_non_constant_vector(self):
        SetHistoricalNonUniformSolutionVector(self.model_part_destination.Nodes, KM.VELOCITY)
        self.mapper.InverseMap(KM.FORCE, KM.VELOCITY)
        mapper_test_case.CheckHistoricalNonUniformValues(self.model_part_origin, KM.FORCE, GetFilePath(self._GetFileName("inverse_map_vector")))
        #mapper_test_case.VtkOutputNodesHistorical(self.model_part_origin, KM.FORCE)

class Projection3D2DMapperNearestElementOrigin2D(mapper_test_case.MapperTestCase):
    '''This class contains the tests for 3D-2D projections for the NearestElement mapper (origin 2D)
    '''

    @classmethod
    def setUpMapper(cls, mapper_parameters):
        mdpa_1 = "3D_blocks_mesh1"
        mdpa_2 = "3D_blocks_mesh2"
        super().setUpModelParts(mdpa_1, mdpa_2)

        # Origin will be 2D
        cls.model_part_origin = cls.model_part_origin.GetSubModelPart("Parts_2D")

        cls.mapper_type = mapper_parameters["mapper_type"].GetString()

        if IsDistributedRun():
            cls.mapper = MappingMPIExtension.MPIMapperFactory.CreateMapper(cls.model_part_origin, cls.model_part_destination, mapper_parameters)
        else:
            cls.mapper = KM.MapperFactory.CreateMapper(cls.model_part_origin, cls.model_part_destination, mapper_parameters)

    @classmethod
    def setUpClass(cls):
        mapper_params = KM.Parameters("""{
            "mapper_type"                        : "projection_3D_2D",
            "base_mapper"                        : "nearest_element",
            "echo_level"                         : 0
        }""")
        cls.setUpMapper(mapper_params)

    def _GetFileName(self, file_appendix):
        return os.path.join("result_files", self.mapper_type, self.__class__.__name__ + "_" + file_appendix)

    def test_Projection3D2DMapper_Map_non_constant_scalar(self):
        SetHistoricalNonUniformSolutionScalar(self.model_part_origin.Nodes, KM.PRESSURE)
        self.mapper.Map(KM.PRESSURE, KM.TEMPERATURE)
        mapper_test_case.CheckHistoricalNonUniformValues(self.model_part_destination, KM.TEMPERATURE, GetFilePath(self._GetFileName("map_scalar")))
        #mapper_test_case.VtkOutputNodesHistorical(self.model_part_destination, KM.TEMPERATURE)

    def test_Projection3D2DMapper_InverseMap_non_constant_scalar(self):
        SetHistoricalNonUniformSolutionScalar(self.model_part_destination.Nodes, KM.TEMPERATURE)
        self.mapper.InverseMap(KM.PRESSURE, KM.TEMPERATURE)
        mapper_test_case.CheckHistoricalNonUniformValues(self.model_part_origin, KM.PRESSURE, GetFilePath(self._GetFileName("inverse_map_scalar")))
        #mapper_test_case.VtkOutputNodesHistorical(self.model_part_origin, KM.PRESSURE)

    def test_Projection3D2DMapper_Map_non_constant_vector(self):
        SetHistoricalNonUniformSolutionVector(self.model_part_origin.Nodes, KM.FORCE)
        self.mapper.Map(KM.FORCE, KM.VELOCITY)
        mapper_test_case.CheckHistoricalNonUniformValues(self.model_part_destination, KM.VELOCITY, GetFilePath(self._GetFileName("map_vector")))
        #mapper_test_case.VtkOutputNodesHistorical(self.model_part_destination, KM.VELOCITY)

    def test_Projection3D2DMapper_InverseMap_non_constant_vector(self):
        SetHistoricalNonUniformSolutionVector(self.model_part_destination.Nodes, KM.VELOCITY)
        self.mapper.InverseMap(KM.FORCE, KM.VELOCITY)
        mapper_test_case.CheckHistoricalNonUniformValues(self.model_part_origin, KM.FORCE, GetFilePath(self._GetFileName("inverse_map_vector")))
        #mapper_test_case.VtkOutputNodesHistorical(self.model_part_origin, KM.FORCE)

class Projection3D2DMapperNearestNeighborDestination2D(mapper_test_case.MapperTestCase):
    '''This class contains the tests for 3D-2D projections for the NearestNeighbor mapper (destination 2D)
    '''

    @classmethod
    def setUpMapper(cls, mapper_parameters):
        mdpa_1 = "3D_blocks_mesh1"
        mdpa_2 = "3D_blocks_mesh2"
        super().setUpModelParts(mdpa_1, mdpa_2)

        # Destination will be 2D
        cls.model_part_destination = cls.model_part_destination.GetSubModelPart("Parts_2D")

        cls.mapper_type = mapper_parameters["mapper_type"].GetString()

        if IsDistributedRun():
            cls.mapper = MappingMPIExtension.MPIMapperFactory.CreateMapper(cls.model_part_origin, cls.model_part_destination, mapper_parameters)
        else:
            cls.mapper = KM.MapperFactory.CreateMapper(cls.model_part_origin, cls.model_part_destination, mapper_parameters)

    @classmethod
    def setUpClass(cls):
        mapper_params = KM.Parameters("""{
            "mapper_type"                        : "projection_3D_2D",
            "base_mapper"                        : "nearest_neighbor",
            "echo_level"                         : 0
        }""")
        cls.setUpMapper(mapper_params)

    def _GetFileName(self, file_appendix):
        return os.path.join("result_files", self.mapper_type, self.__class__.__name__ + "_" + file_appendix)

    def test_Projection3D2DMapper_Map_non_constant_scalar(self):
        SetHistoricalNonUniformSolutionScalar(self.model_part_origin.Nodes, KM.PRESSURE)
        self.mapper.Map(KM.PRESSURE, KM.TEMPERATURE)
        mapper_test_case.CheckHistoricalNonUniformValues(self.model_part_destination, KM.TEMPERATURE, GetFilePath(self._GetFileName("map_scalar")))
        #mapper_test_case.VtkOutputNodesHistorical(self.model_part_destination, KM.TEMPERATURE)

    def test_Projection3D2DMapper_InverseMap_non_constant_scalar(self):
        SetHistoricalNonUniformSolutionScalar(self.model_part_destination.Nodes, KM.TEMPERATURE)
        self.mapper.InverseMap(KM.PRESSURE, KM.TEMPERATURE)
        mapper_test_case.CheckHistoricalNonUniformValues(self.model_part_origin, KM.PRESSURE, GetFilePath(self._GetFileName("inverse_map_scalar")))
        #mapper_test_case.VtkOutputNodesHistorical(self.model_part_origin, KM.PRESSURE)

    def test_Projection3D2DMapper_Map_non_constant_vector(self):
        SetHistoricalNonUniformSolutionVector(self.model_part_origin.Nodes, KM.FORCE)
        self.mapper.Map(KM.FORCE, KM.VELOCITY)
        mapper_test_case.CheckHistoricalNonUniformValues(self.model_part_destination, KM.VELOCITY, GetFilePath(self._GetFileName("map_vector")))
        #mapper_test_case.VtkOutputNodesHistorical(self.model_part_destination, KM.VELOCITY)

    def test_Projection3D2DMapper_InverseMap_non_constant_vector(self):
        SetHistoricalNonUniformSolutionVector(self.model_part_destination.Nodes, KM.VELOCITY)
        self.mapper.InverseMap(KM.FORCE, KM.VELOCITY)
        mapper_test_case.CheckHistoricalNonUniformValues(self.model_part_origin, KM.FORCE, GetFilePath(self._GetFileName("inverse_map_vector")))
        #mapper_test_case.VtkOutputNodesHistorical(self.model_part_origin, KM.FORCE)

class Projection3D2DMapperNearestElementDestination2D(mapper_test_case.MapperTestCase):
    '''This class contains the tests for 3D-2D projections for the NearestElement mapper (destination 2D)
    '''

    @classmethod
    def setUpMapper(cls, mapper_parameters):
        mdpa_1 = "3D_blocks_mesh1"
        mdpa_2 = "3D_blocks_mesh2"
        super().setUpModelParts(mdpa_1, mdpa_2)

        # Destination will be 2D
        cls.model_part_destination = cls.model_part_destination.GetSubModelPart("Parts_2D")

        cls.mapper_type = mapper_parameters["mapper_type"].GetString()

        if IsDistributedRun():
            cls.mapper = MappingMPIExtension.MPIMapperFactory.CreateMapper(cls.model_part_origin, cls.model_part_destination, mapper_parameters)
        else:
            cls.mapper = KM.MapperFactory.CreateMapper(cls.model_part_origin, cls.model_part_destination, mapper_parameters)

    @classmethod
    def setUpClass(cls):
        mapper_params = KM.Parameters("""{
            "mapper_type"                        : "projection_3D_2D",
            "base_mapper"                        : "nearest_element",
            "echo_level"                         : 0
        }""")
        cls.setUpMapper(mapper_params)

    def _GetFileName(self, file_appendix):
        return os.path.join("result_files", self.mapper_type, self.__class__.__name__ + "_" + file_appendix)

    def test_Projection3D2DMapper_Map_non_constant_scalar(self):
        SetHistoricalNonUniformSolutionScalar(self.model_part_origin.Nodes, KM.PRESSURE)
        self.mapper.Map(KM.PRESSURE, KM.TEMPERATURE)
        mapper_test_case.CheckHistoricalNonUniformValues(self.model_part_destination, KM.TEMPERATURE, GetFilePath(self._GetFileName("map_scalar")))
        #mapper_test_case.VtkOutputNodesHistorical(self.model_part_destination, KM.TEMPERATURE)

    def test_Projection3D2DMapper_InverseMap_non_constant_scalar(self):
        SetHistoricalNonUniformSolutionScalar(self.model_part_destination.Nodes, KM.TEMPERATURE)
        self.mapper.InverseMap(KM.PRESSURE, KM.TEMPERATURE)
        mapper_test_case.CheckHistoricalNonUniformValues(self.model_part_origin, KM.PRESSURE, GetFilePath(self._GetFileName("inverse_map_scalar")))
        #mapper_test_case.VtkOutputNodesHistorical(self.model_part_origin, KM.PRESSURE)

    def test_Projection3D2DMapper_Map_non_constant_vector(self):
        SetHistoricalNonUniformSolutionVector(self.model_part_origin.Nodes, KM.FORCE)
        self.mapper.Map(KM.FORCE, KM.VELOCITY)
        mapper_test_case.CheckHistoricalNonUniformValues(self.model_part_destination, KM.VELOCITY, GetFilePath(self._GetFileName("map_vector")))
        #mapper_test_case.VtkOutputNodesHistorical(self.model_part_destination, KM.VELOCITY)

    def test_Projection3D2DMapper_InverseMap_non_constant_vector(self):
        SetHistoricalNonUniformSolutionVector(self.model_part_destination.Nodes, KM.VELOCITY)
        self.mapper.InverseMap(KM.FORCE, KM.VELOCITY)
        mapper_test_case.CheckHistoricalNonUniformValues(self.model_part_origin, KM.FORCE, GetFilePath(self._GetFileName("inverse_map_vector")))
        #mapper_test_case.VtkOutputNodesHistorical(self.model_part_origin, KM.FORCE)

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

if __name__ == '__main__':
    KM.Logger.GetDefaultOutput().SetSeverity(KM.Logger.Severity.WARNING)
    KratosUnittest.main()
