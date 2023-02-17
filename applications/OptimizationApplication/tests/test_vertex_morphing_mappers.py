from typing import Union
from abc import ABC
from abc import abstractmethod

import KratosMultiphysics as Kratos
import KratosMultiphysics.KratosUnittest as kratos_unittest
import KratosMultiphysics.OptimizationApplication as KratosOA
from KratosMultiphysics.OptimizationApplication.utilities.helper_utilities import ContainerVariableDataHolderUnion

class TestVertexMoprhingContainerVariableDataMapper(ABC):
    @classmethod
    def setUpClass(cls):
        cls.KratosEntity = Union[Kratos.Node, Kratos.Condition, Kratos.Element]

        cls.model = Kratos.Model()
        cls.origin_model_part = cls.model.CreateModelPart("origin")
        cls.destination_model_part = cls.model.CreateModelPart("destination")
        cls.origin_model_part.ProcessInfo[Kratos.DOMAIN_SIZE] = 3
        cls.destination_model_part.ProcessInfo[Kratos.DOMAIN_SIZE] = 3

        cls.neighbour_entitys_map = {
            1: [1, 2, 4, 5],
            2: [2, 3, 5, 6],
            3: [4, 5, 7, 8],
            4: [5, 6, 8, 9]
        }

        cls.inverse_neighbour_map = {
            1: [1],
            2: [1, 2],
            3: [2],
            4: [1, 3],
            5: [1, 2, 3, 4],
            6: [2, 4],
            7: [3],
            8: [3, 4],
            9: [4]
        }
        cls.origin_grid_size = 3
        cls.destination_grid_size = 2
        cls.number_of_origin_entities = 9
        cls.number_of_destination_entities = 4

        cls.mapper_parameters = Kratos.Parameters("""{
            "filter_function_type"          : "linear",
            "filter_radius"                 : 1.0,
            "max_entities_in_filter_radius" : 100
        }""")

    @abstractmethod
    def _GetMapper(self) -> Union[KratosOA.Mappers.VertexMorphingNodalContainerVariableDataMapper, KratosOA.Mappers.VertexMorphingConditionContainerVariableDataMapper, KratosOA.Mappers.VertexMorphingElementContainerVariableDataMapper]:
        pass

    @abstractmethod
    def _GetContainerDataHolder(self, model_part: Kratos.ModelPart) -> ContainerVariableDataHolderUnion:
        pass

    @abstractmethod
    def _GetEntity(self, model_part: Kratos.ModelPart, id: int) -> Union[Kratos.Node, Kratos.Condition, Kratos.Element]:
        pass

    @abstractmethod
    def _GetValue(self, entity, variable) -> Union[float, Kratos.Array3]:
        pass

    def test_MapScalar(self):
        mapper = self._GetMapper()
        mapper.Update()

        origin_data = self._GetContainerDataHolder(self.origin_model_part)
        origin_data.ReadDataFromContainerVariable(Kratos.PRESSURE)

        destination_data = self._GetContainerDataHolder(self.destination_model_part)
        mapper.Map(origin_data, destination_data)
        destination_data.AssignDataToContainerVariable(Kratos.PRESSURE)

        for dest_id, orig_ids in self.neighbour_entitys_map.items():
            dest_entity: self.KratosEntity = self._GetEntity(self.destination_model_part, dest_id)
            v = 0.0
            for orig_id in orig_ids:
                orig_entity: self.KratosEntity = self._GetEntity(self.origin_model_part, orig_id)
                v += self._GetValue(orig_entity, Kratos.PRESSURE) * 0.25
            self.assertAlmostEqual(self._GetValue(dest_entity, Kratos.PRESSURE), v, 12)

        inverse_origin_data = self._GetContainerDataHolder(self.origin_model_part)
        mapper.InverseMap(inverse_origin_data, destination_data)
        inverse_origin_data.AssignDataToContainerVariable(Kratos.DENSITY)

        for orig_id, dest_ids in self.inverse_neighbour_map.items():
            orig_entity: self.KratosEntity = self._GetEntity(self.origin_model_part, orig_id)
            v = 0.0
            for dest_id in dest_ids:
                dest_entity: self.KratosEntity = self._GetEntity(self.destination_model_part, dest_id)
                v += self._GetValue(dest_entity, Kratos.PRESSURE) * 0.25
            self.assertAlmostEqual(self._GetValue(orig_entity, Kratos.DENSITY), v, 12)

    def test_MapArray(self):
        mapper = self._GetMapper()
        mapper.Update()

        origin_data = self._GetContainerDataHolder(self.origin_model_part)
        origin_data.ReadDataFromContainerVariable(Kratos.VELOCITY)

        destination_data = self._GetContainerDataHolder(self.destination_model_part)
        mapper.Map(origin_data, destination_data)
        destination_data.AssignDataToContainerVariable(Kratos.VELOCITY)

        for dest_id, orig_ids in self.neighbour_entitys_map.items():
            dest_entity: self.KratosEntity = self._GetEntity(self.destination_model_part, dest_id)
            v = Kratos.Array3([0, 0, 0])
            for orig_id in orig_ids:
                orig_entity: self.KratosEntity = self._GetEntity(self.origin_model_part, orig_id)
                v += self._GetValue(orig_entity, Kratos.VELOCITY) * 0.25
            self.assertVectorAlmostEqual(self._GetValue(dest_entity, Kratos.VELOCITY), v, 12)

        inverse_origin_data = self._GetContainerDataHolder(self.origin_model_part)
        mapper.InverseMap(inverse_origin_data, destination_data)
        inverse_origin_data.AssignDataToContainerVariable(Kratos.ACCELERATION)

        for orig_id, dest_ids in self.inverse_neighbour_map.items():
            orig_entity: self.KratosEntity = self._GetEntity(self.origin_model_part, orig_id)
            v = Kratos.Array3([0, 0, 0])
            for dest_id in dest_ids:
                dest_entity: self.KratosEntity = self._GetEntity(self.destination_model_part, dest_id)
                v += self._GetValue(dest_entity, Kratos.VELOCITY) * 0.25
            self.assertVectorAlmostEqual(self._GetValue(orig_entity, Kratos.ACCELERATION), v, 12)

class TestVertexMorphingNodalContainerVariableDataMapper(TestVertexMoprhingContainerVariableDataMapper, kratos_unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        super().setUpClass()

        cls.origin_model_part.AddNodalSolutionStepVariable(Kratos.PRESSURE)
        cls.origin_model_part.AddNodalSolutionStepVariable(Kratos.DENSITY)
        cls.origin_model_part.AddNodalSolutionStepVariable(Kratos.VELOCITY)
        cls.origin_model_part.AddNodalSolutionStepVariable(Kratos.ACCELERATION)

        cls.destination_model_part.AddNodalSolutionStepVariable(Kratos.PRESSURE)
        cls.destination_model_part.AddNodalSolutionStepVariable(Kratos.DENSITY)
        cls.destination_model_part.AddNodalSolutionStepVariable(Kratos.VELOCITY)
        cls.destination_model_part.AddNodalSolutionStepVariable(Kratos.ACCELERATION)

        for i in range(cls.number_of_origin_entities):
            cls.origin_model_part.CreateNewNode(i+1, i // cls.origin_grid_size, i % cls.origin_grid_size, 0)

        for i in range(cls.number_of_destination_entities):
            cls.destination_model_part.CreateNewNode(i+1, 1.0 / cls.destination_grid_size + (i // cls.destination_grid_size), 1.0 / cls.destination_grid_size + (i % cls.destination_grid_size), 0)

        node: Kratos.Node
        for node in cls.origin_model_part.Nodes:
            node.SetSolutionStepValue(Kratos.PRESSURE, node.Id)
            node.SetSolutionStepValue(Kratos.VELOCITY, Kratos.Array3([node.Id + 1, node.Id + 2, node.Id + 3]))

    def _GetMapper(self):
        return KratosOA.Mappers.VertexMorphingNodalContainerVariableDataMapper(self.origin_model_part, self.destination_model_part, self.mapper_parameters.Clone())

    def _GetContainerDataHolder(self, model_part: Kratos.ModelPart):
        return KratosOA.HistoricalContainerVariableDataHolder(model_part)

    def _GetEntity(self, model_part: Kratos.ModelPart, id: int):
        return model_part.GetNode(id)

    def _GetValue(self, entity, variable):
        return entity.GetSolutionStepValue(variable)

class TestVertexMorphingConditionContainerVariableDataMapper(TestVertexMoprhingContainerVariableDataMapper, kratos_unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        super().setUpClass()

        for i in range(cls.number_of_origin_entities):
            cls.__GenerateEntity(Kratos.Array3([i // cls.origin_grid_size, i % cls.origin_grid_size, 0]), cls.origin_model_part)

        for i in range(cls.number_of_destination_entities):
            cls.__GenerateEntity(Kratos.Array3([1.0 / cls.destination_grid_size + (i // cls.destination_grid_size), 1.0 / cls.destination_grid_size + (i % cls.destination_grid_size), 0]), cls.destination_model_part)

        condition: Kratos.Condition
        for condition in cls.origin_model_part.Conditions:
            condition.SetValue(Kratos.PRESSURE, condition.Id)
            condition.SetValue(Kratos.VELOCITY, Kratos.Array3([condition.Id + 1, condition.Id + 2, condition.Id + 3]))

    def _GetMapper(self):
        return KratosOA.Mappers.VertexMorphingConditionContainerVariableDataMapper(self.origin_model_part, self.destination_model_part, self.mapper_parameters.Clone())

    def _GetContainerDataHolder(self, model_part: Kratos.ModelPart):
        return KratosOA.ConditionContainerVariableDataHolder(model_part)

    def _GetEntity(self, model_part: Kratos.ModelPart, id: int):
        return model_part.GetCondition(id)

    def _GetValue(self, entity, variable):
        return entity.GetValue(variable)

    @classmethod
    def __GenerateEntity(cls, center: Kratos.Array3, model_part: Kratos.ModelPart):
        def __get_index(position):
            for node in model_part.Nodes:
                distance = (position[0] - node.X)**2 + (position[1] - node.Y)**2 + (position[2] - node.Z)**2
                if distance < 1e-8:
                    return node.Id

            node = model_part.CreateNewNode(model_part.NumberOfNodes(), position[0], position[1], position[2])
            return node.Id

        model_part.CreateNewCondition("SurfaceCondition3D4N", model_part.NumberOfConditions() + 1, [
                __get_index(center + Kratos.Array3([-0.5, -0.5, 0.0])),
                __get_index(center + Kratos.Array3([0.5, -0.5, 0.0])),
                __get_index(center + Kratos.Array3([0.5, 0.5, 0.0])),
                __get_index(center + Kratos.Array3([-0.5, 0.5, 0.0]))
        ], model_part.GetProperties()[0])

class TestVertexMorphingElementContainerVariableDataMapper(TestVertexMoprhingContainerVariableDataMapper, kratos_unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        super().setUpClass()

        for i in range(cls.number_of_origin_entities):
            cls.__GenerateEntity(Kratos.Array3([i // cls.origin_grid_size, i % cls.origin_grid_size, 0]), cls.origin_model_part)

        for i in range(cls.number_of_destination_entities):
            cls.__GenerateEntity(Kratos.Array3([1.0 / cls.destination_grid_size + (i // cls.destination_grid_size), 1.0 / cls.destination_grid_size + (i % cls.destination_grid_size), 0]), cls.destination_model_part)

        condition: Kratos.Element
        for condition in cls.origin_model_part.Elements:
            condition.SetValue(Kratos.PRESSURE, condition.Id)
            condition.SetValue(Kratos.VELOCITY, Kratos.Array3([condition.Id + 1, condition.Id + 2, condition.Id + 3]))

    def _GetMapper(self):
        return KratosOA.Mappers.VertexMorphingElementContainerVariableDataMapper(self.origin_model_part, self.destination_model_part, self.mapper_parameters.Clone())

    def _GetContainerDataHolder(self, model_part: Kratos.ModelPart):
        return KratosOA.ElementContainerVariableDataHolder(model_part)

    def _GetEntity(self, model_part: Kratos.ModelPart, id: int):
        return model_part.GetElement(id)

    def _GetValue(self, entity, variable):
        return entity.GetValue(variable)

    @classmethod
    def __GenerateEntity(cls, center: Kratos.Array3, model_part: Kratos.ModelPart):
        def __get_index(position):
            for node in model_part.Nodes:
                distance = (position[0] - node.X)**2 + (position[1] - node.Y)**2 + (position[2] - node.Z)**2
                if distance < 1e-8:
                    return node.Id

            node = model_part.CreateNewNode(model_part.NumberOfNodes(), position[0], position[1], position[2])
            return node.Id

        model_part.CreateNewElement("Element2D4N", model_part.NumberOfElements() + 1, [
                __get_index(center + Kratos.Array3([-0.5, -0.5, 0.0])),
                __get_index(center + Kratos.Array3([0.5, -0.5, 0.0])),
                __get_index(center + Kratos.Array3([0.5, 0.5, 0.0])),
                __get_index(center + Kratos.Array3([-0.5, 0.5, 0.0]))
        ], model_part.GetProperties()[0])

if __name__ == "__main__":
    Kratos.Tester.SetVerbosity(Kratos.Tester.Verbosity.PROGRESS)  # TESTS_OUTPUTS
    kratos_unittest.main()