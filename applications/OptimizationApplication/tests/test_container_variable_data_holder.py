
from abc import ABC
from abc import abstractmethod
import KratosMultiphysics as Kratos
import KratosMultiphysics.OptimizationApplication as KratosOA
from KratosMultiphysics.testing.utilities import ReadModelPart

# Import KratosUnittest
import KratosMultiphysics.KratosUnittest as kratos_unittest

class TestContainerVariableData(ABC):
    @classmethod
    def CreateEntities(cls):
        cls.model = Kratos.Model()
        cls.model_part = cls.model.CreateModelPart("test")
        cls.model_part.AddNodalSolutionStepVariable(Kratos.DENSITY)
        cls.model_part.AddNodalSolutionStepVariable(Kratos.PRESSURE)
        cls.model_part.AddNodalSolutionStepVariable(Kratos.ACCELERATION)
        cls.model_part.AddNodalSolutionStepVariable(Kratos.VELOCITY)
        cls.model_part.ProcessInfo[Kratos.DOMAIN_SIZE] = 3
        ReadModelPart("model_part_utils_test/quads", cls.model_part)

        for node in cls.model_part.Nodes:
            id = node.Id
            node.SetSolutionStepValue(Kratos.VELOCITY, Kratos.Array3([id+3, id+4, id+5]))
            node.SetSolutionStepValue(Kratos.PRESSURE, id+3)
            node.SetValue(Kratos.PRESSURE, id+3)
            node.SetValue(Kratos.VELOCITY, Kratos.Array3([id+3, id+4, id+5]))

        KratosOA.OptimizationUtils.CreateEntitySpecificPropertiesForContainer(cls.model_part, cls.model_part.Conditions)
        for condition in cls.model_part.Conditions:
            id = condition.Id
            condition.Properties[Kratos.PRESSURE] = id+400
            condition.Properties[Kratos.VELOCITY] = Kratos.Array3([id+500, id+600, id+700])
            condition.SetValue(Kratos.PRESSURE, id+4)
            condition.SetValue(Kratos.VELOCITY, Kratos.Array3([id+5, id+6, id+7]))

        KratosOA.OptimizationUtils.CreateEntitySpecificPropertiesForContainer(cls.model_part, cls.model_part.Elements)
        for element in cls.model_part.Elements:
            id = element.Id
            element.Properties[Kratos.PRESSURE] =  id+500
            element.Properties[Kratos.VELOCITY] =  Kratos.Array3([id+600, id+700, id+800])
            element.SetValue(Kratos.PRESSURE, id+5)
            element.SetValue(Kratos.VELOCITY, Kratos.Array3([id+6, id+7, id+8]))

    def test_ContaienrDataAdd(self):
        a = self._GetContainerVariableDataHolder()
        b = self._GetContainerVariableDataHolder()

        a.ReadData(Kratos.VELOCITY)
        b.ReadData(Kratos.VELOCITY)

        c = a + b
        c += b

        c.AssignData(Kratos.ACCELERATION)
        for node in c.GetContainer():
            self.assertVectorAlmostEqual(self._GetValue(node, Kratos.ACCELERATION), self._GetValue(node, Kratos.VELOCITY) * 3, 12)

        c = a + 100.0
        c.AssignData(Kratos.ACCELERATION)
        for node in c.GetContainer():
            self.assertVectorAlmostEqual(self._GetValue(node, Kratos.ACCELERATION), self._GetValue(node, Kratos.VELOCITY) + Kratos.Array3([100.0, 100.0, 100.0]), 12)

        c += 200.0
        c.AssignData(Kratos.ACCELERATION)
        for node in c.GetContainer():
            self.assertVectorAlmostEqual(self._GetValue(node, Kratos.ACCELERATION), self._GetValue(node, Kratos.VELOCITY) + Kratos.Array3([300.0, 300.0, 300.0]), 12)

        a = self._GetContainerVariableDataHolder()
        b = self._GetContainerVariableDataHolder()

        a.ReadData(Kratos.PRESSURE)
        b.ReadData(Kratos.PRESSURE)

        c = a + b
        c += b

        c.AssignData(Kratos.DENSITY)
        for node in c.GetContainer():
            self.assertEqual(self._GetValue(node, Kratos.DENSITY), self._GetValue(node, Kratos.PRESSURE) * 3, 12)

        c = a + 100.0
        c.AssignData(Kratos.DENSITY)
        for node in c.GetContainer():
            self.assertEqual(self._GetValue(node, Kratos.DENSITY), self._GetValue(node, Kratos.PRESSURE) + 100.0, 12)

        c += 100.0
        c.AssignData(Kratos.DENSITY)
        for node in c.GetContainer():
            self.assertEqual(self._GetValue(node, Kratos.DENSITY), self._GetValue(node, Kratos.PRESSURE) + 200.0, 12)

    def test_ContainerVariableDataMultiplyAndSubstract(self):
        a = self._GetContainerVariableDataHolder()
        b = self._GetContainerVariableDataHolder()

        a.ReadData(Kratos.VELOCITY)
        b.ReadData(Kratos.VELOCITY)

        c = a * 4 - b
        c *= 2
        c -= a

        c.AssignData(Kratos.ACCELERATION)
        for node in c.GetContainer():
            self.assertVectorAlmostEqual(self._GetValue(node, Kratos.ACCELERATION), self._GetValue(node, Kratos.VELOCITY) * 5, 12)

        c = a - 100.0
        c.AssignData(Kratos.ACCELERATION)
        for node in c.GetContainer():
            self.assertVectorAlmostEqual(self._GetValue(node, Kratos.ACCELERATION), self._GetValue(node, Kratos.VELOCITY) - Kratos.Array3([100.0, 100.0, 100.0]), 12)

        c -= 200.0
        c.AssignData(Kratos.ACCELERATION)
        for node in c.GetContainer():
            self.assertVectorAlmostEqual(self._GetValue(node, Kratos.ACCELERATION), self._GetValue(node, Kratos.VELOCITY) - Kratos.Array3([300.0, 300.0, 300.0]), 12)

        a = self._GetContainerVariableDataHolder()
        b = self._GetContainerVariableDataHolder()

        a.ReadData(Kratos.PRESSURE)
        b.ReadData(Kratos.PRESSURE)

        c = a * 4 - b
        c *= 2
        c -= a

        c.AssignData(Kratos.DENSITY)
        for node in c.GetContainer():
            self.assertEqual(self._GetValue(node, Kratos.DENSITY), self._GetValue(node, Kratos.PRESSURE) * 5, 12)

        c = a - 100.0
        c.AssignData(Kratos.DENSITY)
        for node in c.GetContainer():
            self.assertEqual(self._GetValue(node, Kratos.DENSITY), self._GetValue(node, Kratos.PRESSURE) - 100.0, 12)

        c -= 100.0
        c.AssignData(Kratos.DENSITY)
        for node in c.GetContainer():
            self.assertEqual(self._GetValue(node, Kratos.DENSITY), self._GetValue(node, Kratos.PRESSURE) - 200.0, 12)

        d = c * a
        d *= b
        d.AssignData(Kratos.DENSITY)
        for node in c.GetContainer():
            self.assertEqual(self._GetValue(node, Kratos.DENSITY), (self._GetValue(node, Kratos.PRESSURE) - 200.0) * self._GetValue(node, Kratos.PRESSURE) ** 2, 12)

        a = self._GetContainerVariableDataHolder()
        a.ReadData(Kratos.VELOCITY)

    def test_ContainerVariableDataDivision(self):
        a = self._GetContainerVariableDataHolder()
        a.ReadData(Kratos.VELOCITY)

        c = a / 2.0
        c /= 2.0

        c.AssignData(Kratos.ACCELERATION)
        for node in c.GetContainer():
            self.assertVectorAlmostEqual(self._GetValue(node, Kratos.ACCELERATION), self._GetValue(node, Kratos.VELOCITY) / 4, 12)

        a = self._GetContainerVariableDataHolder()
        a.ReadData(Kratos.PRESSURE)

        c = a / 2.0
        c /= 2.0

        c.AssignData(Kratos.DENSITY)
        for node in c.GetContainer():
            self.assertEqual(self._GetValue(node, Kratos.DENSITY), self._GetValue(node, Kratos.PRESSURE) / 4, 12)

        d = c / a
        d /= (a * 2)
        d.AssignData(Kratos.DENSITY)
        for node in c.GetContainer():
            self.assertEqual(self._GetValue(node, Kratos.DENSITY), 0.5 * ((self._GetValue(node, Kratos.PRESSURE) / 4) / self._GetValue(node, Kratos.PRESSURE)) / self._GetValue(node, Kratos.PRESSURE) , 12)

        a = self._GetContainerVariableDataHolder()
        a.ReadData(Kratos.VELOCITY)

    def test_ContainerVariableDataPow(self):
        a = self._GetContainerVariableDataHolder()
        a.ReadData(Kratos.VELOCITY)

        c = a ** 2.0
        c **= 2.0

        c.AssignData(Kratos.ACCELERATION)
        for node in c.GetContainer():
            ref_value = self._GetValue(node, Kratos.VELOCITY)
            self.assertVectorAlmostEqual(self._GetValue(node, Kratos.ACCELERATION), Kratos.Array3([ref_value[0]**4, ref_value[1]**4, ref_value[2]**4]), 12)

        a = self._GetContainerVariableDataHolder()
        a.ReadData(Kratos.PRESSURE)

        c = a ** 2.0
        c **= 2.0

        c.AssignData(Kratos.DENSITY)
        for node in c.GetContainer():
            self.assertEqual(self._GetValue(node, Kratos.DENSITY), self._GetValue(node, Kratos.PRESSURE) ** 4, 12)

    def test_ContainerVariableDataNeg(self):
        a = self._GetContainerVariableDataHolder()
        a.ReadData(Kratos.VELOCITY)

        c = -a

        c.AssignData(Kratos.ACCELERATION)
        for node in c.GetContainer():
            self.assertVectorAlmostEqual(self._GetValue(node, Kratos.ACCELERATION), self._GetValue(node, Kratos.VELOCITY) * (-1.0), 12)

        a = self._GetContainerVariableDataHolder()
        a.ReadData(Kratos.PRESSURE)

        c = -a

        c.AssignData(Kratos.DENSITY)
        for node in c.GetContainer():
            self.assertEqual(self._GetValue(node, Kratos.DENSITY), self._GetValue(node, Kratos.PRESSURE) * (-1.0), 12)

    def test_SetDataForContainerVariable(self):
        a = self._GetContainerVariableDataHolder()
        a.SetData(Kratos.Array3([1, 2, 3]))
        a.AssignData(Kratos.ACCELERATION)
        for node in a.GetContainer():
            self.assertVectorAlmostEqual(self._GetValue(node, Kratos.ACCELERATION), Kratos.Array3([1, 2, 3]), 12)

        a = self._GetContainerVariableDataHolder()
        a.SetData(10)
        a.AssignData(Kratos.DENSITY)
        for node in a.GetContainer():
            self.assertEqual(self._GetValue(node, Kratos.DENSITY), 10)

    def test_Clone(self):
        a = self._GetContainerVariableDataHolder()

        a.ReadData(Kratos.VELOCITY)

        b = a.Clone()
        b.SetData(Kratos.Array3([10, 11, 12]))

        a.AssignData(Kratos.ACCELERATION)
        for node in a.GetContainer():
            self.assertVectorAlmostEqual(self._GetValue(node, Kratos.ACCELERATION), self._GetValue(node, Kratos.VELOCITY), 12)

        b.AssignData(Kratos.ACCELERATION)
        for node in a.GetContainer():
            self.assertVectorAlmostEqual(self._GetValue(node, Kratos.ACCELERATION), Kratos.Array3([10, 11, 12]), 12)

        a = self._GetContainerVariableDataHolder()

        a.ReadData(Kratos.PRESSURE)

        b = a.Clone()
        b.SetData(12)

        a.AssignData(Kratos.DENSITY)
        for node in a.GetContainer():
            self.assertEqual(self._GetValue(node, Kratos.DENSITY), self._GetValue(node, Kratos.PRESSURE), 12)

        b.AssignData(Kratos.DENSITY)
        for node in a.GetContainer():
            self.assertEqual(self._GetValue(node, Kratos.DENSITY), 12, 12)

    def test_GetContainer(self):
        a = self._GetContainerVariableDataHolder()
        self.assertEqual(self._GetContainer(), a.GetContainer())

    @abstractmethod
    def _GetContainerVariableDataHolder(self):
        pass

    @abstractmethod
    def _GetContainer(self):
        pass

    @abstractmethod
    def _GetValue(self, entity, variable):
        pass

class TestHistoricalContainerVariableDataHolder(kratos_unittest.TestCase, TestContainerVariableData):
    @classmethod
    def setUpClass(cls):
        cls.CreateEntities()

    def test_CopyData(self):
        a = self._GetContainerVariableDataHolder()
        b = KratosOA.NodalNonHistoricalVariableData(self.model_part)

        a.ReadData(Kratos.VELOCITY)
        b.CopyDataFrom(a)

        b.AssignData(Kratos.ACCELERATION)
        for node in b.GetContainer():
            self.assertVectorAlmostEqual(node.GetValue(Kratos.ACCELERATION), self._GetValue(node, Kratos.VELOCITY), 12)

        a = self._GetContainerVariableDataHolder()

        a.ReadData(Kratos.PRESSURE)
        b.CopyDataFrom(a)

        b.AssignData(Kratos.DENSITY)
        for node in b.GetContainer():
            self.assertEqual(node.GetValue(Kratos.DENSITY), self._GetValue(node, Kratos.PRESSURE), 12)

        b = KratosOA.NodalNonHistoricalVariableData(a)
        b += 1
        b.AssignData(Kratos.DENSITY)
        for node in b.GetContainer():
            self.assertEqual(node.GetValue(Kratos.DENSITY), self._GetValue(node, Kratos.PRESSURE) + 1, 12)

    def _GetContainerVariableDataHolder(self):
        return KratosOA.HistoricalVariableData(self.model_part)

    def _GetContainer(self):
        return self.model_part.GetCommunicator().LocalMesh().Nodes

    def _GetValue(self, entity, variable):
        return entity.GetSolutionStepValue(variable)

class TestNodalContainerVariableDataHolder(kratos_unittest.TestCase, TestContainerVariableData):
    @classmethod
    def setUpClass(cls):
        cls.CreateEntities()

    def test_CopyData(self):
        a = self._GetContainerVariableDataHolder()
        b = KratosOA.HistoricalVariableData(self.model_part)

        a.ReadData(Kratos.VELOCITY)
        b.CopyDataFrom(a)

        b.AssignData(Kratos.ACCELERATION)
        for node in b.GetContainer():
            self.assertVectorAlmostEqual(node.GetSolutionStepValue(Kratos.ACCELERATION), self._GetValue(node, Kratos.VELOCITY), 12)

        a = self._GetContainerVariableDataHolder()

        a.ReadData(Kratos.PRESSURE)
        b.CopyDataFrom(a)

        b.AssignData(Kratos.DENSITY)
        for node in b.GetContainer():
            self.assertEqual(node.GetSolutionStepValue(Kratos.DENSITY), self._GetValue(node, Kratos.PRESSURE), 12)

        b = KratosOA.HistoricalVariableData(a)
        b += 1
        b.AssignData(Kratos.DENSITY)
        for node in b.GetContainer():
            self.assertEqual(node.GetSolutionStepValue(Kratos.DENSITY), self._GetValue(node, Kratos.PRESSURE) + 1, 12)

    def _GetContainerVariableDataHolder(self):
        return KratosOA.NodalNonHistoricalVariableData(self.model_part)

    def _GetContainer(self):
        return self.model_part.GetCommunicator().LocalMesh().Nodes

    def _GetValue(self, entity, variable):
        return entity.GetValue(variable)

class TestConditionContainerVariableDataHolder(kratos_unittest.TestCase, TestContainerVariableData):
    @classmethod
    def setUpClass(cls):
        cls.CreateEntities()

    def test_CopyData(self):
        a = self._GetContainerVariableDataHolder()
        b = KratosOA.ConditionPropertiesVariableData(self.model_part)

        a.ReadData(Kratos.VELOCITY)
        b.CopyDataFrom(a)

        b.AssignData(Kratos.ACCELERATION)
        for node in b.GetContainer():
            self.assertVectorAlmostEqual(node.Properties[Kratos.ACCELERATION], self._GetValue(node, Kratos.VELOCITY), 12)

        a = self._GetContainerVariableDataHolder()

        a.ReadData(Kratos.PRESSURE)
        b.CopyDataFrom(a)

        b.AssignData(Kratos.DENSITY)
        for node in b.GetContainer():
            self.assertEqual(node.Properties[Kratos.DENSITY], self._GetValue(node, Kratos.PRESSURE), 12)

        b = KratosOA.ConditionPropertiesVariableData(a)
        b += 1
        b.AssignData(Kratos.DENSITY)
        for node in b.GetContainer():
            self.assertEqual(node.Properties[Kratos.DENSITY], self._GetValue(node, Kratos.PRESSURE) + 1, 12)

    def _GetContainerVariableDataHolder(self):
        return KratosOA.ConditionNonHistoricalVariableData(self.model_part)

    def _GetContainer(self):
        return self.model_part.GetCommunicator().LocalMesh().Conditions

    def _GetValue(self, entity, variable):
        return entity.GetValue(variable)

class TestElementContainerVariableDataHolder(kratos_unittest.TestCase, TestContainerVariableData):
    @classmethod
    def setUpClass(cls):
        cls.CreateEntities()

    def test_CopyData(self):
        a = self._GetContainerVariableDataHolder()
        b = KratosOA.ElementPropertiesVariableData(self.model_part)

        a.ReadData(Kratos.VELOCITY)
        b.CopyDataFrom(a)

        b.AssignData(Kratos.ACCELERATION)
        for node in b.GetContainer():
            self.assertVectorAlmostEqual(node.Properties[Kratos.ACCELERATION], self._GetValue(node, Kratos.VELOCITY), 12)

        a = self._GetContainerVariableDataHolder()

        a.ReadData(Kratos.PRESSURE)
        b.CopyDataFrom(a)

        b.AssignData(Kratos.DENSITY)
        for node in b.GetContainer():
            self.assertEqual(node.Properties[Kratos.DENSITY], self._GetValue(node, Kratos.PRESSURE), 12)

        b = KratosOA.ElementPropertiesVariableData(a)
        b += 1
        b.AssignData(Kratos.DENSITY)
        for node in b.GetContainer():
            self.assertEqual(node.Properties[Kratos.DENSITY], self._GetValue(node, Kratos.PRESSURE) + 1, 12)

    def _GetContainerVariableDataHolder(self):
        return KratosOA.ElementNonHistoricalVariableData(self.model_part)

    def _GetContainer(self):
        return self.model_part.GetCommunicator().LocalMesh().Elements

    def _GetValue(self, entity, variable):
        return entity.GetValue(variable)

class TestConditionPropertiesContainerVariableDataHolder(kratos_unittest.TestCase, TestContainerVariableData):
    @classmethod
    def setUpClass(cls):
        cls.CreateEntities()

    def test_CopyData(self):
        a = self._GetContainerVariableDataHolder()
        b = KratosOA.ConditionNonHistoricalVariableData(self.model_part)

        a.ReadData(Kratos.VELOCITY)
        b.CopyDataFrom(a)

        b.AssignData(Kratos.ACCELERATION)
        for node in b.GetContainer():
            self.assertVectorAlmostEqual(node.GetValue(Kratos.ACCELERATION), self._GetValue(node, Kratos.VELOCITY), 12)

        a = self._GetContainerVariableDataHolder()

        a.ReadData(Kratos.PRESSURE)
        b.CopyDataFrom(a)

        b.AssignData(Kratos.DENSITY)
        for node in b.GetContainer():
            self.assertEqual(node.GetValue(Kratos.DENSITY), self._GetValue(node, Kratos.PRESSURE), 12)

        b = KratosOA.ConditionNonHistoricalVariableData(a)
        b += 1
        b.AssignData(Kratos.DENSITY)
        for node in b.GetContainer():
            self.assertEqual(node.GetValue(Kratos.DENSITY), self._GetValue(node, Kratos.PRESSURE) + 1, 12)

    def _GetContainerVariableDataHolder(self):
        return KratosOA.ConditionPropertiesVariableData(self.model_part)

    def _GetContainer(self):
        return self.model_part.GetCommunicator().LocalMesh().Conditions

    def _GetValue(self, entity, variable):
        return entity.Properties[variable]

class TestElementPropertiesContainerVariableDataHolder(kratos_unittest.TestCase, TestContainerVariableData):
    @classmethod
    def setUpClass(cls):
        cls.CreateEntities()

    def test_CopyData(self):
        a = self._GetContainerVariableDataHolder()
        b = KratosOA.ElementNonHistoricalVariableData(self.model_part)

        a.ReadData(Kratos.VELOCITY)
        b.CopyDataFrom(a)

        b.AssignData(Kratos.ACCELERATION)
        for node in b.GetContainer():
            self.assertVectorAlmostEqual(node.GetValue(Kratos.ACCELERATION), self._GetValue(node, Kratos.VELOCITY), 12)

        a = self._GetContainerVariableDataHolder()

        a.ReadData(Kratos.PRESSURE)
        b.CopyDataFrom(a)

        b.AssignData(Kratos.DENSITY)
        for node in b.GetContainer():
            self.assertEqual(node.GetValue(Kratos.DENSITY), self._GetValue(node, Kratos.PRESSURE), 12)

        b = KratosOA.ElementNonHistoricalVariableData(a)
        b += 1
        b.AssignData(Kratos.DENSITY)
        for node in b.GetContainer():
            self.assertEqual(node.GetValue(Kratos.DENSITY), self._GetValue(node, Kratos.PRESSURE) + 1, 12)

    def _GetContainerVariableDataHolder(self):
        return KratosOA.ElementPropertiesVariableData(self.model_part)

    def _GetContainer(self):
        return self.model_part.GetCommunicator().LocalMesh().Elements

    def _GetValue(self, entity, variable):
        return entity.Properties[variable]

if __name__ == "__main__":
    Kratos.Tester.SetVerbosity(Kratos.Tester.Verbosity.PROGRESS)  # TESTS_OUTPUTS
    kratos_unittest.main()