import KratosMultiphysics as KM
import KratosMultiphysics.MappingApplication as KratosMapping
data_comm = KM.Testing.GetDefaultDataCommunicator()
if data_comm.IsDistributed():
    from KratosMultiphysics.MappingApplication import MPIExtension as MappingMPIExtension

import mapper_test_case
from math import sin, cos
import os

def GetFilePath(file_name):
    return os.path.join(os.path.dirname(os.path.realpath(__file__)), file_name)

class BasicMapperTests(mapper_test_case.MapperTestCase):
    '''This class contains basic tests that every mapper should pass
    This included e.g. testing if mapping a constant field works
    Also it is checked if the mapper-flags are working correctly
    '''

    @classmethod
    def setUpMapper(cls, mapper_parameters, switch_sides=False):
        if switch_sides:
            super().setUpModelParts("cube_quad", "cube_tri")
        else:
            super().setUpModelParts("cube_tri", "cube_quad")
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


    def test_Map_constant_scalar(self):
        val = 1.234
        KM.VariableUtils().SetScalarVar(KM.PRESSURE, val, self.interface_model_part_origin.Nodes)
        self.mapper.Map(KM.PRESSURE, KM.TEMPERATURE)
        self._CheckHistoricalUniformValuesScalar(GetNodes(self.interface_model_part_destination), KM.TEMPERATURE, val)

    def test_InverseMap_constant_scalar(self):
        val = -571.147
        KM.VariableUtils().SetScalarVar(KM.TEMPERATURE, val, self.interface_model_part_destination.Nodes)
        self.mapper.InverseMap(KM.PRESSURE, KM.TEMPERATURE)
        self._CheckHistoricalUniformValuesScalar(GetNodes(self.interface_model_part_origin), KM.PRESSURE, val)

    def test_Map_constant_vector(self):
        val = KM.Vector([1.234, -22.845, 11.775])
        KM.VariableUtils().SetVectorVar(KM.FORCE, val, self.interface_model_part_origin.Nodes)
        self.mapper.Map(KM.FORCE, KM.VELOCITY)
        self._CheckHistoricalUniformValuesVector(GetNodes(self.interface_model_part_destination), KM.VELOCITY, val)

    def test_InverseMap_constant_vector(self):
        val = KM.Vector([-51.234, -22.845, 118.775])
        KM.VariableUtils().SetVectorVar(KM.VELOCITY, val, self.interface_model_part_destination.Nodes)
        self.mapper.InverseMap(KM.FORCE, KM.VELOCITY)
        self._CheckHistoricalUniformValuesVector(GetNodes(self.interface_model_part_origin), KM.FORCE, val)

    def test_Map_non_constant_scalar(self):
        SetHistoricalNonUniformSolutionScalar(self.interface_model_part_origin.Nodes, KM.PRESSURE)
        self.mapper.Map(KM.PRESSURE, KM.TEMPERATURE)
        mapper_test_case.CheckHistoricalNonUniformValues(self.interface_model_part_destination, KM.TEMPERATURE, GetFilePath(self._GetFileName("map_scalar")))

    def test_InverseMap_non_constant_scalar(self):
        SetHistoricalNonUniformSolutionScalar(self.interface_model_part_destination.Nodes, KM.TEMPERATURE)
        self.mapper.InverseMap(KM.PRESSURE, KM.TEMPERATURE)
        mapper_test_case.CheckHistoricalNonUniformValues(self.interface_model_part_origin, KM.PRESSURE, GetFilePath(self._GetFileName("inverse_map_scalar")))

    def test_Map_non_constant_vector(self):
        SetHistoricalNonUniformSolutionVector(self.interface_model_part_origin.Nodes, KM.FORCE)
        self.mapper.Map(KM.FORCE, KM.VELOCITY)
        mapper_test_case.CheckHistoricalNonUniformValues(self.interface_model_part_destination, KM.VELOCITY, GetFilePath(self._GetFileName("map_vector")))

    def test_InverseMap_non_constant_vector(self):
        SetHistoricalNonUniformSolutionVector(self.interface_model_part_destination.Nodes, KM.VELOCITY)
        self.mapper.InverseMap(KM.FORCE, KM.VELOCITY)
        mapper_test_case.CheckHistoricalNonUniformValues(self.interface_model_part_origin, KM.FORCE, GetFilePath(self._GetFileName("inverse_map_vector")))

    def test_SWAP_SIGN_Map_scalar(self):
        val = 1.234
        KM.VariableUtils().SetScalarVar(KM.PRESSURE, val, self.interface_model_part_origin.Nodes)
        self.mapper.Map(KM.PRESSURE, KM.TEMPERATURE, KM.Mapper.SWAP_SIGN)
        self._CheckHistoricalUniformValuesScalar(GetNodes(self.interface_model_part_destination), KM.TEMPERATURE, -val)

    def test_SWAP_SIGN_InverseMap_scalar(self):
        val = -571.147
        KM.VariableUtils().SetScalarVar(KM.TEMPERATURE, val, self.interface_model_part_destination.Nodes)
        self.mapper.InverseMap(KM.PRESSURE, KM.TEMPERATURE, KM.Mapper.SWAP_SIGN)
        self._CheckHistoricalUniformValuesScalar(GetNodes(self.interface_model_part_origin), KM.PRESSURE, -val)

    def test_SWAP_SIGN_Map_vector(self):
        val = KM.Vector([1.234, -22.845, 11.775])
        KM.VariableUtils().SetVectorVar(KM.FORCE, val, self.interface_model_part_origin.Nodes)
        self.mapper.Map(KM.FORCE, KM.VELOCITY, KM.Mapper.SWAP_SIGN)
        self._CheckHistoricalUniformValuesVector(GetNodes(self.interface_model_part_destination), KM.VELOCITY, [(-1)*x for x in val])

    def test_SWAP_SIGN_InverseMap_vector(self):
        val = KM.Vector([-51.234, -22.845, 118.775])
        KM.VariableUtils().SetVectorVar(KM.VELOCITY, val, self.interface_model_part_destination.Nodes)
        self.mapper.InverseMap(KM.FORCE, KM.VELOCITY, KM.Mapper.SWAP_SIGN)
        self._CheckHistoricalUniformValuesVector(GetNodes(self.interface_model_part_origin), KM.FORCE, [(-1)*x for x in val])

    def test_ADD_VALUES_Map_scalar(self):
        val_1 = 1.234
        val_2 = -571.147
        KM.VariableUtils().SetScalarVar(KM.PRESSURE, val_1, self.interface_model_part_origin.Nodes)
        self.mapper.Map(KM.PRESSURE, KM.TEMPERATURE) # set the initial field

        KM.VariableUtils().SetScalarVar(KM.PRESSURE, val_2, self.interface_model_part_origin.Nodes)
        self.mapper.Map(KM.PRESSURE, KM.TEMPERATURE, KM.Mapper.ADD_VALUES)
        self._CheckHistoricalUniformValuesScalar(GetNodes(self.interface_model_part_destination), KM.TEMPERATURE, val_1+val_2)

    def test_ADD_VALUES_InverseMap_scalar(self):
        val_1 = -571.147
        val_2 = 128.336
        KM.VariableUtils().SetScalarVar(KM.TEMPERATURE, val_1, self.interface_model_part_destination.Nodes)
        self.mapper.InverseMap(KM.PRESSURE, KM.TEMPERATURE)

        KM.VariableUtils().SetScalarVar(KM.TEMPERATURE, val_2, self.interface_model_part_destination.Nodes)
        self.mapper.InverseMap(KM.PRESSURE, KM.TEMPERATURE, KM.Mapper.ADD_VALUES)
        self._CheckHistoricalUniformValuesScalar(GetNodes(self.interface_model_part_origin), KM.PRESSURE, val_1+val_2)

    def test_ADD_VALUES_Map_vector(self):
        val_1 = KM.Vector([1.234, -22.845, 11.83])
        val_2 = KM.Vector([-51.9234, -22.845, 118.775])
        KM.VariableUtils().SetVectorVar(KM.FORCE, val_1, self.interface_model_part_origin.Nodes)
        self.mapper.Map(KM.FORCE, KM.VELOCITY) # set the initial field

        KM.VariableUtils().SetVectorVar(KM.FORCE, val_2, self.interface_model_part_origin.Nodes)
        self.mapper.Map(KM.FORCE, KM.VELOCITY, KM.Mapper.ADD_VALUES)
        self._CheckHistoricalUniformValuesVector(GetNodes(self.interface_model_part_destination), KM.VELOCITY, val_1+val_2)

    def test_ADD_VALUES_InverseMap_vector(self):
        val_1 = KM.Vector([1.234, -22.845, 11.83])
        val_2 = KM.Vector([-51.9234, -22.845, 118.775])
        KM.VariableUtils().SetVectorVar(KM.VELOCITY, val_1, self.interface_model_part_destination.Nodes)
        self.mapper.InverseMap(KM.FORCE, KM.VELOCITY) # set the initial field

        KM.VariableUtils().SetVectorVar(KM.VELOCITY, val_2, self.interface_model_part_destination.Nodes)
        self.mapper.InverseMap(KM.FORCE, KM.VELOCITY, KM.Mapper.ADD_VALUES)
        self._CheckHistoricalUniformValuesVector(GetNodes(self.interface_model_part_origin), KM.FORCE, val_1+val_2)

    def test_SWAP_SIGN_and_ADD_VALUES_scalar(self):
        val_1 = 1.234
        val_2 = -571.147
        KM.VariableUtils().SetScalarVar(KM.PRESSURE, val_1, self.interface_model_part_origin.Nodes)
        self.mapper.Map(KM.PRESSURE, KM.TEMPERATURE) # set the initial field

        KM.VariableUtils().SetScalarVar(KM.PRESSURE, val_2, self.interface_model_part_origin.Nodes)
        self.mapper.Map(KM.PRESSURE, KM.TEMPERATURE, KM.Mapper.ADD_VALUES | KM.Mapper.SWAP_SIGN)
        self._CheckHistoricalUniformValuesScalar(GetNodes(self.interface_model_part_destination), KM.TEMPERATURE, val_1-val_2)

    def test_Map_USE_TRANSPOSE_constant_scalar(self):
        val = 1.234
        KM.VariableUtils().SetScalarVar(KM.PRESSURE, val, self.interface_model_part_origin.Nodes)
        self.mapper.Map(KM.PRESSURE, KM.TEMPERATURE, KM.Mapper.USE_TRANSPOSE)

        sum_origin = KM.VariableUtils().SumHistoricalNodeScalarVariable(KM.PRESSURE, self.interface_model_part_origin, 0)
        sum_destination = KM.VariableUtils().SumHistoricalNodeScalarVariable(KM.TEMPERATURE, self.interface_model_part_destination, 0)
        self.assertAlmostEqual(sum_origin, sum_destination)

    def test_InverseMap_USE_TRANSPOSE_constant_scalar(self):
        val = 1.234
        KM.VariableUtils().SetScalarVar(KM.TEMPERATURE, val, self.interface_model_part_destination.Nodes)
        self.mapper.InverseMap(KM.PRESSURE, KM.TEMPERATURE, KM.Mapper.USE_TRANSPOSE)

        sum_origin = KM.VariableUtils().SumHistoricalNodeScalarVariable(KM.PRESSURE, self.interface_model_part_origin, 0)
        sum_destination = KM.VariableUtils().SumHistoricalNodeScalarVariable(KM.TEMPERATURE, self.interface_model_part_destination, 0)
        self.assertAlmostEqual(sum_origin, sum_destination)

    def test_Map_USE_TRANSPOSE_constant_vector(self):
        val = KM.Vector([1.234, -22.845, 11.83])
        KM.VariableUtils().SetVectorVar(KM.FORCE, val, self.interface_model_part_origin.Nodes)
        self.mapper.Map(KM.FORCE, KM.VELOCITY, KM.Mapper.USE_TRANSPOSE)

        sum_origin = KM.VariableUtils().SumHistoricalNodeVectorVariable(KM.FORCE, self.interface_model_part_origin, 0)
        sum_destination = KM.VariableUtils().SumHistoricalNodeVectorVariable(KM.VELOCITY, self.interface_model_part_destination, 0)
        self.assertVectorAlmostEqual(sum_origin, sum_destination)

    def test_InverseMap_USE_TRANSPOSE_constant_vector(self):
        val = KM.Vector([1.234, -22.845, 11.83])
        KM.VariableUtils().SetVectorVar(KM.VELOCITY, val, self.interface_model_part_destination.Nodes)
        self.mapper.InverseMap(KM.FORCE, KM.VELOCITY, KM.Mapper.USE_TRANSPOSE)

        sum_origin = KM.VariableUtils().SumHistoricalNodeVectorVariable(KM.FORCE, self.interface_model_part_origin, 0)
        sum_destination = KM.VariableUtils().SumHistoricalNodeVectorVariable(KM.VELOCITY, self.interface_model_part_destination, 0)
        self.assertVectorAlmostEqual(sum_origin, sum_destination)

    def test_Map_constant_scalar_TO_NON_HISTORICAL(self):
        val = 9.234
        KM.VariableUtils().SetVariable(KM.PRESSURE, val, self.interface_model_part_origin.Nodes)
        self.mapper.Map(KM.PRESSURE, KM.TEMPERATURE, KM.Mapper.TO_NON_HISTORICAL)
        self._CheckUniformValuesScalar(GetNodes(self.interface_model_part_destination), KM.TEMPERATURE, val)

    def test_Map_constant_scalar_FROM_NON_HISTORICAL(self):
        val = -961.234
        KM.VariableUtils().SetNonHistoricalVariable(KM.PRESSURE, val, self.interface_model_part_origin.Nodes)
        self.mapper.Map(KM.PRESSURE, KM.TEMPERATURE, KM.Mapper.FROM_NON_HISTORICAL)
        self._CheckHistoricalUniformValuesScalar(GetNodes(self.interface_model_part_destination), KM.TEMPERATURE, val)

    def test_Map_constant_scalar_both_non_historical(self):
        val = 34.234
        KM.VariableUtils().SetNonHistoricalVariable(KM.PRESSURE, val, self.interface_model_part_origin.Nodes)
        self.mapper.Map(KM.PRESSURE, KM.TEMPERATURE, KM.Mapper.FROM_NON_HISTORICAL | KM.Mapper.TO_NON_HISTORICAL)
        self._CheckUniformValuesScalar(GetNodes(self.interface_model_part_destination), KM.TEMPERATURE, val)

    def test_InverseMap_constant_scalar_TO_NON_HISTORICAL(self):
        val = 8.23554
        KM.VariableUtils().SetVariable(KM.TEMPERATURE, val, self.interface_model_part_destination.Nodes)
        self.mapper.InverseMap(KM.PRESSURE, KM.TEMPERATURE, KM.Mapper.TO_NON_HISTORICAL)
        self._CheckUniformValuesScalar(GetNodes(self.interface_model_part_origin), KM.PRESSURE, val)

    def test_InverseMap_constant_scalar_FROM_NON_HISTORICAL(self):
        val = -96741.234
        KM.VariableUtils().SetNonHistoricalVariable(KM.TEMPERATURE, val, self.interface_model_part_destination.Nodes)
        self.mapper.InverseMap(KM.PRESSURE, KM.TEMPERATURE, KM.Mapper.FROM_NON_HISTORICAL)
        self._CheckHistoricalUniformValuesScalar(GetNodes(self.interface_model_part_origin), KM.PRESSURE, val)

    def test_InverseMap_constant_scalar_both_non_historical(self):
        val = 3134.24734
        KM.VariableUtils().SetNonHistoricalVariable(KM.TEMPERATURE, val, self.interface_model_part_destination.Nodes)
        self.mapper.InverseMap(KM.PRESSURE, KM.TEMPERATURE, KM.Mapper.FROM_NON_HISTORICAL | KM.Mapper.TO_NON_HISTORICAL)
        self._CheckUniformValuesScalar(GetNodes(self.interface_model_part_origin), KM.PRESSURE, val)


    def test_Map_USE_TRANSPOSE_constant_scalar_TO_NON_HISTORICAL(self):
        val = 17.09
        KM.VariableUtils().SetScalarVar(KM.PRESSURE, val, self.interface_model_part_origin.Nodes)
        self.mapper.Map(KM.PRESSURE, KM.TEMPERATURE, KM.Mapper.USE_TRANSPOSE | KM.Mapper.TO_NON_HISTORICAL)

        sum_origin = KM.VariableUtils().SumHistoricalNodeScalarVariable(KM.PRESSURE, self.interface_model_part_origin, 0)
        sum_destination = KM.VariableUtils().SumNonHistoricalNodeScalarVariable(KM.TEMPERATURE, self.interface_model_part_destination)
        self.assertAlmostEqual(sum_origin, sum_destination)

    def test_Map_USE_TRANSPOSE_constant_scalar_FROM_NON_HISTORICAL(self):
        val = -88.76
        KM.VariableUtils().SetNonHistoricalVariable(KM.PRESSURE, val, self.interface_model_part_origin.Nodes)
        self.mapper.Map(KM.PRESSURE, KM.TEMPERATURE, KM.Mapper.USE_TRANSPOSE | KM.Mapper.FROM_NON_HISTORICAL)

        sum_origin = KM.VariableUtils().SumNonHistoricalNodeScalarVariable(KM.PRESSURE, self.interface_model_part_origin)
        sum_destination = KM.VariableUtils().SumHistoricalNodeScalarVariable(KM.TEMPERATURE, self.interface_model_part_destination, 0)
        self.assertAlmostEqual(sum_origin, sum_destination)

    def test_Map_USE_TRANSPOSE_constant_scalar_both_non_historical(self):
        val = 101.234
        KM.VariableUtils().SetNonHistoricalVariable(KM.PRESSURE, val, self.interface_model_part_origin.Nodes)
        self.mapper.Map(KM.PRESSURE, KM.TEMPERATURE, KM.Mapper.USE_TRANSPOSE | KM.Mapper.FROM_NON_HISTORICAL | KM.Mapper.TO_NON_HISTORICAL)

        sum_origin = KM.VariableUtils().SumNonHistoricalNodeScalarVariable(KM.PRESSURE, self.interface_model_part_origin)
        sum_destination = KM.VariableUtils().SumNonHistoricalNodeScalarVariable(KM.TEMPERATURE, self.interface_model_part_destination)
        self.assertAlmostEqual(sum_origin, sum_destination)

    def test_InverseMap_USE_TRANSPOSE_constant_scalar_TO_NON_HISTORICAL(self):
        val = 23.189
        KM.VariableUtils().SetScalarVar(KM.TEMPERATURE, val, self.interface_model_part_destination.Nodes)
        self.mapper.InverseMap(KM.PRESSURE, KM.TEMPERATURE, KM.Mapper.USE_TRANSPOSE | KM.Mapper.TO_NON_HISTORICAL)

        sum_origin = KM.VariableUtils().SumNonHistoricalNodeScalarVariable(KM.PRESSURE, self.interface_model_part_origin)
        sum_destination = KM.VariableUtils().SumHistoricalNodeScalarVariable(KM.TEMPERATURE, self.interface_model_part_destination, 0)
        self.assertAlmostEqual(sum_origin, sum_destination)

    def test_InverseMap_USE_TRANSPOSE_constant_scalar_FROM_NON_HISTORICAL(self):
        val = 651.234
        KM.VariableUtils().SetNonHistoricalVariable(KM.TEMPERATURE, val, self.interface_model_part_destination.Nodes)
        self.mapper.InverseMap(KM.PRESSURE, KM.TEMPERATURE, KM.Mapper.USE_TRANSPOSE | KM.Mapper.FROM_NON_HISTORICAL)

        sum_origin = KM.VariableUtils().SumHistoricalNodeScalarVariable(KM.PRESSURE, self.interface_model_part_origin, 0)
        sum_destination = KM.VariableUtils().SumNonHistoricalNodeScalarVariable(KM.TEMPERATURE, self.interface_model_part_destination)
        self.assertAlmostEqual(sum_origin, sum_destination)

    def test_InverseMap_USE_TRANSPOSE_constant_scalar_both_non_historical(self):
        val = 53.761
        KM.VariableUtils().SetNonHistoricalVariable(KM.TEMPERATURE, val, self.interface_model_part_destination.Nodes)
        self.mapper.InverseMap(KM.PRESSURE, KM.TEMPERATURE, KM.Mapper.USE_TRANSPOSE | KM.Mapper.FROM_NON_HISTORICAL | KM.Mapper.TO_NON_HISTORICAL)

        sum_origin = KM.VariableUtils().SumNonHistoricalNodeScalarVariable(KM.PRESSURE, self.interface_model_part_origin)
        sum_destination = KM.VariableUtils().SumNonHistoricalNodeScalarVariable(KM.TEMPERATURE, self.interface_model_part_destination)
        self.assertAlmostEqual(sum_origin, sum_destination)


    # def test_UpdateInterface(self):
    #     pass

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
        val = 12*sin(node.X0) + node.Y0*15 + 22*node.Z0
        node.SetSolutionStepValue(variable, val)

def SetHistoricalNonUniformSolutionVector(nodes, variable):
    for node in nodes:
        val_1 = 12*sin(node.X0) + node.Y0*15 + 22*node.Z0
        val_2 = 33*cos(node.X0) + node.Y0*5 + 22*node.Z0
        val_3 = 12*sin(node.Y0) + node.Z0*15 + 22*node.X0
        node.SetSolutionStepValue(variable, KM.Vector([val_1, val_2, val_3]))

def GetNodes(model_part):
    return model_part.Nodes
