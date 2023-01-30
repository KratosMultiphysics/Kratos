
from abc import ABC
from abc import abstractmethod
import KratosMultiphysics as Kratos
import KratosMultiphysics.OptimizationApplication as KratosOA

# Import KratosUnittest
import KratosMultiphysics.KratosUnittest as kratos_unittest

class TestContainerDataBase(ABC):
    @classmethod
    def CreateEntities(cls):
        cls.model = Kratos.Model()
        cls.model_part = cls.model.CreateModelPart("test")
        cls.model_part.AddNodalSolutionStepVariable(Kratos.DENSITY)
        cls.model_part.AddNodalSolutionStepVariable(Kratos.PRESSURE)
        cls.model_part.AddNodalSolutionStepVariable(Kratos.ACCELERATION)
        cls.model_part.AddNodalSolutionStepVariable(Kratos.VELOCITY)
        cls.model_part.ProcessInfo[Kratos.DOMAIN_SIZE] = 3

        number_of_nodes = 10
        for id in range(1, number_of_nodes + 1):
            node = cls.model_part.CreateNewNode(id, id, id+1, id+2)
            node.SetSolutionStepValue(Kratos.VELOCITY, Kratos.Array3([id+3, id+4, id+5]))
            node.SetSolutionStepValue(Kratos.PRESSURE, id+3)
            node.SetValue(Kratos.PRESSURE, id+3)
            node.SetValue(Kratos.VELOCITY, Kratos.Array3([id+3, id+4, id+5]))

        number_of_conditions = 11
        for id in range(1, number_of_conditions + 1):
            properties = cls.model_part.CreateNewProperties(id)
            properties.SetValue(Kratos.PRESSURE, id+400)
            properties.SetValue(Kratos.VELOCITY, Kratos.Array3([id+500, id+600, id+700]))
            condition = cls.model_part.CreateNewCondition("LineCondition2D2N", id + 1, [(id % number_of_nodes) + 1, ((id + 1) % number_of_nodes) + 1 ], properties)
            condition.SetValue(Kratos.PRESSURE, id+4)
            condition.SetValue(Kratos.VELOCITY, Kratos.Array3([id+5, id+6, id+7]))

        number_of_elements = 12
        for id in range(1, number_of_elements + 1):
            properties = cls.model_part.CreateNewProperties(id + number_of_conditions)
            properties.SetValue(Kratos.PRESSURE, id+500)
            properties.SetValue(Kratos.VELOCITY, Kratos.Array3([id+600, id+700, id+800]))
            element = cls.model_part.CreateNewElement("Element2D3N", id + 2, [(id % number_of_nodes) + 1, ((id + 1) % number_of_nodes) + 1, ((id + 2) % number_of_nodes) + 1 ], properties)
            element.SetValue(Kratos.PRESSURE, id+5)
            element.SetValue(Kratos.VELOCITY, Kratos.Array3([id+6, id+7, id+8]))

    def test_ContaienrDataAdd(self):
        a = self._GetContainerData()
        b = self._GetContainerData()

        a.ReadDataFromContainerVariable(Kratos.VELOCITY)
        b.ReadDataFromContainerVariable(Kratos.VELOCITY)

        c = a + b
        c += b

        c.AssignDataToContainerVariable(Kratos.ACCELERATION)
        for node in c.GetContainer():
            self.assertVectorAlmostEqual(self._GetValue(node, Kratos.ACCELERATION), self._GetValue(node, Kratos.VELOCITY) * 3, 12)

        c = a + 100.0
        c.AssignDataToContainerVariable(Kratos.ACCELERATION)
        for node in c.GetContainer():
            self.assertVectorAlmostEqual(self._GetValue(node, Kratos.ACCELERATION), self._GetValue(node, Kratos.VELOCITY) + Kratos.Array3([100.0, 100.0, 100.0]), 12)

        c += 200.0
        c.AssignDataToContainerVariable(Kratos.ACCELERATION)
        for node in c.GetContainer():
            self.assertVectorAlmostEqual(self._GetValue(node, Kratos.ACCELERATION), self._GetValue(node, Kratos.VELOCITY) + Kratos.Array3([300.0, 300.0, 300.0]), 12)

        a = self._GetContainerData()
        b = self._GetContainerData()

        a.ReadDataFromContainerVariable(Kratos.PRESSURE)
        b.ReadDataFromContainerVariable(Kratos.PRESSURE)

        c = a + b
        c += b

        c.AssignDataToContainerVariable(Kratos.DENSITY)
        for node in c.GetContainer():
            self.assertEqual(self._GetValue(node, Kratos.DENSITY), self._GetValue(node, Kratos.PRESSURE) * 3, 12)

        c = a + 100.0
        c.AssignDataToContainerVariable(Kratos.DENSITY)
        for node in c.GetContainer():
            self.assertEqual(self._GetValue(node, Kratos.DENSITY), self._GetValue(node, Kratos.PRESSURE) + 100.0, 12)

        c += 100.0
        c.AssignDataToContainerVariable(Kratos.DENSITY)
        for node in c.GetContainer():
            self.assertEqual(self._GetValue(node, Kratos.DENSITY), self._GetValue(node, Kratos.PRESSURE) + 200.0, 12)

    def test_ContainerDataMultiplyAndSubstract(self):
        a = self._GetContainerData()
        b = self._GetContainerData()

        a.ReadDataFromContainerVariable(Kratos.VELOCITY)
        b.ReadDataFromContainerVariable(Kratos.VELOCITY)

        c = a * 4 - b
        c *= 2
        c -= a

        c.AssignDataToContainerVariable(Kratos.ACCELERATION)
        for node in c.GetContainer():
            self.assertVectorAlmostEqual(self._GetValue(node, Kratos.ACCELERATION), self._GetValue(node, Kratos.VELOCITY) * 5, 12)

        c = a - 100.0
        c.AssignDataToContainerVariable(Kratos.ACCELERATION)
        for node in c.GetContainer():
            self.assertVectorAlmostEqual(self._GetValue(node, Kratos.ACCELERATION), self._GetValue(node, Kratos.VELOCITY) - Kratos.Array3([100.0, 100.0, 100.0]), 12)

        c -= 200.0
        c.AssignDataToContainerVariable(Kratos.ACCELERATION)
        for node in c.GetContainer():
            self.assertVectorAlmostEqual(self._GetValue(node, Kratos.ACCELERATION), self._GetValue(node, Kratos.VELOCITY) - Kratos.Array3([300.0, 300.0, 300.0]), 12)

        a = self._GetContainerData()
        b = self._GetContainerData()

        a.ReadDataFromContainerVariable(Kratos.PRESSURE)
        b.ReadDataFromContainerVariable(Kratos.PRESSURE)

        c = a * 4 - b
        c *= 2
        c -= a

        c.AssignDataToContainerVariable(Kratos.DENSITY)
        for node in c.GetContainer():
            self.assertEqual(self._GetValue(node, Kratos.DENSITY), self._GetValue(node, Kratos.PRESSURE) * 5, 12)

        c = a - 100.0
        c.AssignDataToContainerVariable(Kratos.DENSITY)
        for node in c.GetContainer():
            self.assertEqual(self._GetValue(node, Kratos.DENSITY), self._GetValue(node, Kratos.PRESSURE) - 100.0, 12)

        c -= 100.0
        c.AssignDataToContainerVariable(Kratos.DENSITY)
        for node in c.GetContainer():
            self.assertEqual(self._GetValue(node, Kratos.DENSITY), self._GetValue(node, Kratos.PRESSURE) - 200.0, 12)

    def test_ContainerDataDivision(self):
        a = self._GetContainerData()
        a.ReadDataFromContainerVariable(Kratos.VELOCITY)

        c = a / 2.0
        c /= 2.0

        c.AssignDataToContainerVariable(Kratos.ACCELERATION)
        for node in c.GetContainer():
            self.assertVectorAlmostEqual(self._GetValue(node, Kratos.ACCELERATION), self._GetValue(node, Kratos.VELOCITY) / 4, 12)

        a = self._GetContainerData()
        a.ReadDataFromContainerVariable(Kratos.PRESSURE)

        c = a / 2.0
        c /= 2.0

        c.AssignDataToContainerVariable(Kratos.DENSITY)
        for node in c.GetContainer():
            self.assertEqual(self._GetValue(node, Kratos.DENSITY), self._GetValue(node, Kratos.PRESSURE) / 4, 12)

    def test_ContainerDataPow(self):
        a = self._GetContainerData()
        a.ReadDataFromContainerVariable(Kratos.VELOCITY)

        c = a ** 2.0
        c **= 2.0

        c.AssignDataToContainerVariable(Kratos.ACCELERATION)
        for node in c.GetContainer():
            ref_value = self._GetValue(node, Kratos.VELOCITY)
            self.assertVectorAlmostEqual(self._GetValue(node, Kratos.ACCELERATION), Kratos.Array3([ref_value[0]**4, ref_value[1]**4, ref_value[2]**4]), 12)

        a = self._GetContainerData()
        a.ReadDataFromContainerVariable(Kratos.PRESSURE)

        c = a ** 2.0
        c **= 2.0

        c.AssignDataToContainerVariable(Kratos.DENSITY)
        for node in c.GetContainer():
            self.assertEqual(self._GetValue(node, Kratos.DENSITY), self._GetValue(node, Kratos.PRESSURE) ** 4, 12)

    def test_ContainerDataNeg(self):
        a = self._GetContainerData()
        a.ReadDataFromContainerVariable(Kratos.VELOCITY)

        c = -a

        c.AssignDataToContainerVariable(Kratos.ACCELERATION)
        for node in c.GetContainer():
            self.assertVectorAlmostEqual(self._GetValue(node, Kratos.ACCELERATION), self._GetValue(node, Kratos.VELOCITY) * (-1.0), 12)

        a = self._GetContainerData()
        a.ReadDataFromContainerVariable(Kratos.PRESSURE)

        c = -a

        c.AssignDataToContainerVariable(Kratos.DENSITY)
        for node in c.GetContainer():
            self.assertEqual(self._GetValue(node, Kratos.DENSITY), self._GetValue(node, Kratos.PRESSURE) * (-1.0), 12)

    def test_SetDataForContainerVariable(self):
        a = self._GetContainerData()
        a.SetDataForContainerVariable(Kratos.VELOCITY, Kratos.Array3([1, 2, 3]))
        a.AssignDataToContainerVariable(Kratos.ACCELERATION)
        for node in a.GetContainer():
            self.assertVectorAlmostEqual(self._GetValue(node, Kratos.ACCELERATION), Kratos.Array3([1, 2, 3]), 12)

        a = self._GetContainerData()
        a.SetDataForContainerVariable(Kratos.PRESSURE, 10)
        a.AssignDataToContainerVariable(Kratos.DENSITY)
        for node in a.GetContainer():
            self.assertEqual(self._GetValue(node, Kratos.DENSITY), 10)

    def test_Clone(self):
        a = self._GetContainerData()

        a.ReadDataFromContainerVariable(Kratos.VELOCITY)

        b = a.Clone()
        b.SetDataForContainerVariable(Kratos.VELOCITY, Kratos.Array3([10, 11, 12]))

        a.AssignDataToContainerVariable(Kratos.ACCELERATION)
        for node in a.GetContainer():
            self.assertVectorAlmostEqual(self._GetValue(node, Kratos.ACCELERATION), self._GetValue(node, Kratos.VELOCITY), 12)

        b.AssignDataToContainerVariable(Kratos.ACCELERATION)
        for node in a.GetContainer():
            self.assertVectorAlmostEqual(self._GetValue(node, Kratos.ACCELERATION), Kratos.Array3([10, 11, 12]), 12)

        a = self._GetContainerData()

        a.ReadDataFromContainerVariable(Kratos.PRESSURE)

        b = a.Clone()
        b.SetDataForContainerVariable(Kratos.PRESSURE, 12)

        a.AssignDataToContainerVariable(Kratos.DENSITY)
        for node in a.GetContainer():
            self.assertEqual(self._GetValue(node, Kratos.DENSITY), self._GetValue(node, Kratos.PRESSURE), 12)

        b.AssignDataToContainerVariable(Kratos.DENSITY)
        for node in a.GetContainer():
            self.assertEqual(self._GetValue(node, Kratos.DENSITY), 12, 12)

    def test_GetContainer(self):
        a = self._GetContainerData()
        self.assertEqual(self._GetContainer(), a.GetContainer())

    @abstractmethod
    def _GetContainerData(self):
        pass

    @abstractmethod
    def _GetContainer(self):
        pass

    @abstractmethod
    def _GetValue(self, entity, variable):
        pass

class TestHistoricalContainerData(kratos_unittest.TestCase, TestContainerDataBase):
    @classmethod
    def setUpClass(cls):
        cls.CreateEntities()

    def test_CopyData(self):
        a = self._GetContainerData()
        b = KratosOA.NodalContainerData(self.model_part)

        a.ReadDataFromContainerVariable(Kratos.VELOCITY)
        b.CopyDataFrom(a)

        b.AssignDataToContainerVariable(Kratos.ACCELERATION)
        for node in b.GetContainer():
            self.assertVectorAlmostEqual(node.GetValue(Kratos.ACCELERATION), self._GetValue(node, Kratos.VELOCITY), 12)

        a = self._GetContainerData()
        b = KratosOA.NodalContainerData(self.model_part)

        a.ReadDataFromContainerVariable(Kratos.PRESSURE)
        b.CopyDataFrom(a)

        b.AssignDataToContainerVariable(Kratos.DENSITY)
        for node in b.GetContainer():
            self.assertEqual(node.GetValue(Kratos.DENSITY), self._GetValue(node, Kratos.PRESSURE), 12)

    def _GetContainerData(self):
        return KratosOA.HistoricalContainerData(self.model_part)

    def _GetContainer(self):
        return self.model_part.Nodes

    def _GetValue(self, entity, variable):
        return entity.GetSolutionStepValue(variable)

class TestNodalContainerData(kratos_unittest.TestCase, TestContainerDataBase):
    @classmethod
    def setUpClass(cls):
        cls.CreateEntities()

    def test_CopyData(self):
        a = self._GetContainerData()
        b = KratosOA.HistoricalContainerData(self.model_part)

        a.ReadDataFromContainerVariable(Kratos.VELOCITY)
        b.CopyDataFrom(a)

        b.AssignDataToContainerVariable(Kratos.ACCELERATION)
        for node in b.GetContainer():
            self.assertVectorAlmostEqual(node.GetSolutionStepValue(Kratos.ACCELERATION), self._GetValue(node, Kratos.VELOCITY), 12)

        a = self._GetContainerData()
        b = KratosOA.HistoricalContainerData(self.model_part)

        a.ReadDataFromContainerVariable(Kratos.PRESSURE)
        b.CopyDataFrom(a)

        b.AssignDataToContainerVariable(Kratos.DENSITY)
        for node in b.GetContainer():
            self.assertEqual(node.GetSolutionStepValue(Kratos.DENSITY), self._GetValue(node, Kratos.PRESSURE), 12)

    def _GetContainerData(self):
        return KratosOA.NodalContainerData(self.model_part)

    def _GetContainer(self):
        return self.model_part.Nodes

    def _GetValue(self, entity, variable):
        return entity.GetValue(variable)

class TestConditionContainerData(kratos_unittest.TestCase, TestContainerDataBase):
    @classmethod
    def setUpClass(cls):
        cls.CreateEntities()

    def test_CopyData(self):
        a = self._GetContainerData()
        b = KratosOA.ConditionPropertiesContainerData(self.model_part)

        a.ReadDataFromContainerVariable(Kratos.VELOCITY)
        b.CopyDataFrom(a)

        b.AssignDataToContainerVariable(Kratos.ACCELERATION)
        for node in b.GetContainer():
            self.assertVectorAlmostEqual(node.Properties[Kratos.ACCELERATION], self._GetValue(node, Kratos.VELOCITY), 12)

        a = self._GetContainerData()
        b = KratosOA.ConditionPropertiesContainerData(self.model_part)

        a.ReadDataFromContainerVariable(Kratos.PRESSURE)
        b.CopyDataFrom(a)

        b.AssignDataToContainerVariable(Kratos.DENSITY)
        for node in b.GetContainer():
            self.assertEqual(node.Properties[Kratos.DENSITY], self._GetValue(node, Kratos.PRESSURE), 12)

    def _GetContainerData(self):
        return KratosOA.ConditionContainerData(self.model_part)

    def _GetContainer(self):
        return self.model_part.Conditions

    def _GetValue(self, entity, variable):
        return entity.GetValue(variable)

class TestElementContainerData(kratos_unittest.TestCase, TestContainerDataBase):
    @classmethod
    def setUpClass(cls):
        cls.CreateEntities()

    def test_CopyData(self):
        a = self._GetContainerData()
        b = KratosOA.ElementPropertiesContainerData(self.model_part)

        a.ReadDataFromContainerVariable(Kratos.VELOCITY)
        b.CopyDataFrom(a)

        b.AssignDataToContainerVariable(Kratos.ACCELERATION)
        for node in b.GetContainer():
            self.assertVectorAlmostEqual(node.Properties[Kratos.ACCELERATION], self._GetValue(node, Kratos.VELOCITY), 12)

        a = self._GetContainerData()
        b = KratosOA.ElementPropertiesContainerData(self.model_part)

        a.ReadDataFromContainerVariable(Kratos.PRESSURE)
        b.CopyDataFrom(a)

        b.AssignDataToContainerVariable(Kratos.DENSITY)
        for node in b.GetContainer():
            self.assertEqual(node.Properties[Kratos.DENSITY], self._GetValue(node, Kratos.PRESSURE), 12)

    def _GetContainerData(self):
        return KratosOA.ElementContainerData(self.model_part)

    def _GetContainer(self):
        return self.model_part.Elements

    def _GetValue(self, entity, variable):
        return entity.GetValue(variable)

class TestConditionPropertiesContainerData(kratos_unittest.TestCase, TestContainerDataBase):
    @classmethod
    def setUpClass(cls):
        cls.CreateEntities()

    def test_CopyData(self):
        a = self._GetContainerData()
        b = KratosOA.ConditionContainerData(self.model_part)

        a.ReadDataFromContainerVariable(Kratos.VELOCITY)
        b.CopyDataFrom(a)

        b.AssignDataToContainerVariable(Kratos.ACCELERATION)
        for node in b.GetContainer():
            self.assertVectorAlmostEqual(node.GetValue(Kratos.ACCELERATION), self._GetValue(node, Kratos.VELOCITY), 12)

        a = self._GetContainerData()
        b = KratosOA.ConditionContainerData(self.model_part)

        a.ReadDataFromContainerVariable(Kratos.PRESSURE)
        b.CopyDataFrom(a)

        b.AssignDataToContainerVariable(Kratos.DENSITY)
        for node in b.GetContainer():
            self.assertEqual(node.GetValue(Kratos.DENSITY), self._GetValue(node, Kratos.PRESSURE), 12)

    def _GetContainerData(self):
        return KratosOA.ConditionPropertiesContainerData(self.model_part)

    def _GetContainer(self):
        return self.model_part.Conditions

    def _GetValue(self, entity, variable):
        return entity.Properties[variable]

class TestElementPropertiesContainerData(kratos_unittest.TestCase, TestContainerDataBase):
    @classmethod
    def setUpClass(cls):
        cls.CreateEntities()

    def test_CopyData(self):
        a = self._GetContainerData()
        b = KratosOA.ElementContainerData(self.model_part)

        a.ReadDataFromContainerVariable(Kratos.VELOCITY)
        b.CopyDataFrom(a)

        b.AssignDataToContainerVariable(Kratos.ACCELERATION)
        for node in b.GetContainer():
            self.assertVectorAlmostEqual(node.GetValue(Kratos.ACCELERATION), self._GetValue(node, Kratos.VELOCITY), 12)

        a = self._GetContainerData()
        b = KratosOA.ElementContainerData(self.model_part)

        a.ReadDataFromContainerVariable(Kratos.PRESSURE)
        b.CopyDataFrom(a)

        b.AssignDataToContainerVariable(Kratos.DENSITY)
        for node in b.GetContainer():
            self.assertEqual(node.GetValue(Kratos.DENSITY), self._GetValue(node, Kratos.PRESSURE), 12)

    def _GetContainerData(self):
        return KratosOA.ElementPropertiesContainerData(self.model_part)

    def _GetContainer(self):
        return self.model_part.Elements

    def _GetValue(self, entity, variable):
        return entity.Properties[variable]

if __name__ == "__main__":
    Kratos.Tester.SetVerbosity(Kratos.Tester.Verbosity.PROGRESS)  # TESTS_OUTPUTS
    kratos_unittest.main()