from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

import KratosMultiphysics as KM
import KratosMultiphysics.MappingApplication as KratosMapping
data_comm = KM.DataCommunicator.GetDefault()
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
            super(BasicMapperTests, cls).setUpModelParts("cube_quad", "cube_tri")
        else:
            super(BasicMapperTests, cls).setUpModelParts("cube_tri", "cube_quad")
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
            cls.mapper = KratosMapping.MapperFactory.CreateMPIMapper(
                cls.model_part_origin, cls.model_part_destination, mapper_parameters)
        else:
            cls.mapper = KratosMapping.MapperFactory.CreateMapper(
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
        self.mapper.Map(KM.PRESSURE, KM.TEMPERATURE, KratosMapping.Mapper.SWAP_SIGN)
        self._CheckHistoricalUniformValuesScalar(GetNodes(self.interface_model_part_destination), KM.TEMPERATURE, -val)

    def test_SWAP_SIGN_InverseMap_scalar(self):
        val = -571.147
        KM.VariableUtils().SetScalarVar(KM.TEMPERATURE, val, self.interface_model_part_destination.Nodes)
        self.mapper.InverseMap(KM.PRESSURE, KM.TEMPERATURE, KratosMapping.Mapper.SWAP_SIGN)
        self._CheckHistoricalUniformValuesScalar(GetNodes(self.interface_model_part_origin), KM.PRESSURE, -val)

    def test_SWAP_SIGN_Map_vector(self):
        val = KM.Vector([1.234, -22.845, 11.775])
        KM.VariableUtils().SetVectorVar(KM.FORCE, val, self.interface_model_part_origin.Nodes)
        self.mapper.Map(KM.FORCE, KM.VELOCITY, KratosMapping.Mapper.SWAP_SIGN)
        self._CheckHistoricalUniformValuesVector(GetNodes(self.interface_model_part_destination), KM.VELOCITY, [(-1)*x for x in val])

    def test_SWAP_SIGN_InverseMap_vector(self):
        val = KM.Vector([-51.234, -22.845, 118.775])
        KM.VariableUtils().SetVectorVar(KM.VELOCITY, val, self.interface_model_part_destination.Nodes)
        self.mapper.InverseMap(KM.FORCE, KM.VELOCITY, KratosMapping.Mapper.SWAP_SIGN)
        self._CheckHistoricalUniformValuesVector(GetNodes(self.interface_model_part_origin), KM.FORCE, [(-1)*x for x in val])

    def test_ADD_VALUES_Map_scalar(self):
        val_1 = 1.234
        val_2 = -571.147
        KM.VariableUtils().SetScalarVar(KM.PRESSURE, val_1, self.interface_model_part_origin.Nodes)
        self.mapper.Map(KM.PRESSURE, KM.TEMPERATURE) # set the initial field

        KM.VariableUtils().SetScalarVar(KM.PRESSURE, val_2, self.interface_model_part_origin.Nodes)
        self.mapper.Map(KM.PRESSURE, KM.TEMPERATURE, KratosMapping.Mapper.ADD_VALUES)
        self._CheckHistoricalUniformValuesScalar(GetNodes(self.interface_model_part_destination), KM.TEMPERATURE, val_1+val_2)

    def test_ADD_VALUES_InverseMap_scalar(self):
        val_1 = -571.147
        val_2 = 128.336
        KM.VariableUtils().SetScalarVar(KM.TEMPERATURE, val_1, self.interface_model_part_destination.Nodes)
        self.mapper.InverseMap(KM.PRESSURE, KM.TEMPERATURE)

        KM.VariableUtils().SetScalarVar(KM.TEMPERATURE, val_2, self.interface_model_part_destination.Nodes)
        self.mapper.InverseMap(KM.PRESSURE, KM.TEMPERATURE, KratosMapping.Mapper.ADD_VALUES)
        self._CheckHistoricalUniformValuesScalar(GetNodes(self.interface_model_part_origin), KM.PRESSURE, val_1+val_2)

    def test_ADD_VALUES_Map_vector(self):
        val_1 = KM.Vector([1.234, -22.845, 11.83])
        val_2 = KM.Vector([-51.9234, -22.845, 118.775])
        KM.VariableUtils().SetVectorVar(KM.FORCE, val_1, self.interface_model_part_origin.Nodes)
        self.mapper.Map(KM.FORCE, KM.VELOCITY) # set the initial field

        KM.VariableUtils().SetVectorVar(KM.FORCE, val_2, self.interface_model_part_origin.Nodes)
        self.mapper.Map(KM.FORCE, KM.VELOCITY, KratosMapping.Mapper.ADD_VALUES)
        self._CheckHistoricalUniformValuesVector(GetNodes(self.interface_model_part_destination), KM.VELOCITY, val_1+val_2)

    def test_ADD_VALUES_InverseMap_vector(self):
        val_1 = KM.Vector([1.234, -22.845, 11.83])
        val_2 = KM.Vector([-51.9234, -22.845, 118.775])
        KM.VariableUtils().SetVectorVar(KM.VELOCITY, val_1, self.interface_model_part_destination.Nodes)
        self.mapper.InverseMap(KM.FORCE, KM.VELOCITY) # set the initial field

        KM.VariableUtils().SetVectorVar(KM.VELOCITY, val_2, self.interface_model_part_destination.Nodes)
        self.mapper.InverseMap(KM.FORCE, KM.VELOCITY, KratosMapping.Mapper.ADD_VALUES)
        self._CheckHistoricalUniformValuesVector(GetNodes(self.interface_model_part_origin), KM.FORCE, val_1+val_2)

    def test_SWAP_SIGN_and_ADD_VALUES_scalar(self):
        val_1 = 1.234
        val_2 = -571.147
        KM.VariableUtils().SetScalarVar(KM.PRESSURE, val_1, self.interface_model_part_origin.Nodes)
        self.mapper.Map(KM.PRESSURE, KM.TEMPERATURE) # set the initial field

        KM.VariableUtils().SetScalarVar(KM.PRESSURE, val_2, self.interface_model_part_origin.Nodes)
        self.mapper.Map(KM.PRESSURE, KM.TEMPERATURE, KratosMapping.Mapper.ADD_VALUES | KratosMapping.Mapper.SWAP_SIGN)
        self._CheckHistoricalUniformValuesScalar(GetNodes(self.interface_model_part_destination), KM.TEMPERATURE, val_1-val_2)

    def test_Map_USE_TRANSPOSE_constant_scalar(self):
        val = 1.234
        KM.VariableUtils().SetScalarVar(KM.PRESSURE, val, self.interface_model_part_origin.Nodes)
        self.mapper.Map(KM.PRESSURE, KM.TEMPERATURE, KratosMapping.Mapper.USE_TRANSPOSE)

        sum_origin = KM.VariableUtils().SumHistoricalNodeScalarVariable(KM.PRESSURE, self.interface_model_part_origin, 0)
        sum_destination = KM.VariableUtils().SumHistoricalNodeScalarVariable(KM.TEMPERATURE, self.interface_model_part_destination, 0)
        self.assertAlmostEqual(sum_origin, sum_destination)

    def test_InverseMap_USE_TRANSPOSE_constant_scalar(self):
        val = 1.234
        KM.VariableUtils().SetScalarVar(KM.TEMPERATURE, val, self.interface_model_part_destination.Nodes)
        self.mapper.InverseMap(KM.PRESSURE, KM.TEMPERATURE, KratosMapping.Mapper.USE_TRANSPOSE)

        sum_origin = KM.VariableUtils().SumHistoricalNodeScalarVariable(KM.PRESSURE, self.interface_model_part_origin, 0)
        sum_destination = KM.VariableUtils().SumHistoricalNodeScalarVariable(KM.TEMPERATURE, self.interface_model_part_destination, 0)
        self.assertAlmostEqual(sum_origin, sum_destination)

    def test_Map_USE_TRANSPOSE_constant_vector(self):
        val = KM.Vector([1.234, -22.845, 11.83])
        KM.VariableUtils().SetVectorVar(KM.FORCE, val, self.interface_model_part_origin.Nodes)
        self.mapper.Map(KM.FORCE, KM.VELOCITY, KratosMapping.Mapper.USE_TRANSPOSE)

        sum_origin = KM.VariableUtils().SumHistoricalNodeVectorVariable(KM.FORCE, self.interface_model_part_origin, 0)
        sum_destination = KM.VariableUtils().SumHistoricalNodeVectorVariable(KM.VELOCITY, self.interface_model_part_destination, 0)
        self.assertAlmostEqual(sum_origin[0], sum_destination[0])
        self.assertAlmostEqual(sum_origin[1], sum_destination[1])
        self.assertAlmostEqual(sum_origin[2], sum_destination[2])

    def test_InverseMap_USE_TRANSPOSE_constant_vector(self):
        val = KM.Vector([1.234, -22.845, 11.83])
        KM.VariableUtils().SetVectorVar(KM.VELOCITY, val, self.interface_model_part_destination.Nodes)
        self.mapper.InverseMap(KM.FORCE, KM.VELOCITY, KratosMapping.Mapper.USE_TRANSPOSE)

        sum_origin = KM.VariableUtils().SumHistoricalNodeVectorVariable(KM.FORCE, self.interface_model_part_origin, 0)
        sum_destination = KM.VariableUtils().SumHistoricalNodeVectorVariable(KM.VELOCITY, self.interface_model_part_destination, 0)
        self.assertAlmostEqual(sum_origin[0], sum_destination[0])
        self.assertAlmostEqual(sum_origin[1], sum_destination[1])
        self.assertAlmostEqual(sum_origin[2], sum_destination[2])

    def test_Is_conforming(self):
        is_conforming = self.mapper.AreMeshesConforming()
        self.assertEqual(is_conforming, True)

    def test_Is_not_conforming(self):
        non_conform_parameters = self.mapper_parameters.Clone()
        non_conform_parameters.AddEmptyValue("search_radius").SetDouble(1e-6)

        if data_comm.IsDistributed():
            map_creator = KratosMapping.MapperFactory.CreateMPIMapper
        else:
            map_creator = KratosMapping.MapperFactory.CreateMapper

        non_conform_mapper = map_creator(
            self.model_part_origin,
            self.model_part_destination,
            non_conform_parameters
        )

        is_conforming = non_conform_mapper.AreMeshesConforming()
        self.assertEqual(is_conforming, False)


    # def test_UpdateInterface(self):
    #     pass

    # def test_TO_NON_HISTORICAL(self):
    #     pass

    # def test_FROM_NON_HISTORICAL(self):
    #     pass

    # def test_both_NON_HISTORICAL(self):
    #     pass



    def _CheckHistoricalUniformValuesScalar(self, nodes, variable, exp_value):
        for node in nodes:
            self.assertAlmostEqual(node.GetSolutionStepValue(variable), exp_value)

    def _CheckHistoricalUniformValuesVector(self, nodes, variable, exp_value):
        for node in nodes:
            nodal_val = node.GetSolutionStepValue(variable)
            self.assertAlmostEqual(nodal_val[0], exp_value[0])
            self.assertAlmostEqual(nodal_val[1], exp_value[1])
            self.assertAlmostEqual(nodal_val[2], exp_value[2])

    def _CheckUniformValuesScalar(self, entities, variable, exp_value):
        for entity in entities:
            self.assertAlmostEqual(entity.GetValue(variable), exp_value)

    def _CheckUniformValuesVector(self, entities, variable, exp_value):
        for entity in entities:
            val = entity.GetValue(variable)
            self.assertAlmostEqual(val[0], exp_value[0])
            self.assertAlmostEqual(val[1], exp_value[1])
            self.assertAlmostEqual(val[2], exp_value[2])

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
    return model_part.GetCommunicator().LocalMesh().Nodes
    # return model_part.Nodes # TODO this is the correct version, requires some synchronization though!
