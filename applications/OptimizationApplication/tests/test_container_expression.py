
from abc import ABC
from abc import abstractmethod
import KratosMultiphysics as Kratos
import KratosMultiphysics.OptimizationApplication as KratosOA
from KratosMultiphysics.testing.utilities import ReadModelPart

# Import KratosUnittest
import KratosMultiphysics.KratosUnittest as kratos_unittest

class TestContainerExpression(ABC):
    @classmethod
    def CreateEntities(cls):
        cls.model = Kratos.Model()
        cls.model_part = cls.model.CreateModelPart("test")
        cls.model_part.AddNodalSolutionStepVariable(Kratos.DENSITY)
        cls.model_part.AddNodalSolutionStepVariable(Kratos.PRESSURE)
        cls.model_part.AddNodalSolutionStepVariable(Kratos.ACCELERATION)
        cls.model_part.AddNodalSolutionStepVariable(Kratos.VELOCITY)
        with kratos_unittest.WorkFolderScope(".", __file__, True):
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
        a = self._GetContainerExpression()
        b = self._GetContainerExpression()

        a.Read(Kratos.VELOCITY)
        b.Read(Kratos.VELOCITY)

        c = a + b
        c += b

        c.Evaluate(Kratos.ACCELERATION)
        for node in c.GetContainer():
            self.assertVectorAlmostEqual(self._GetValue(node, Kratos.ACCELERATION), self._GetValue(node, Kratos.VELOCITY) * 3, 12)

        c = a + 100.0
        c.Evaluate(Kratos.ACCELERATION)
        for node in c.GetContainer():
            self.assertVectorAlmostEqual(self._GetValue(node, Kratos.ACCELERATION), self._GetValue(node, Kratos.VELOCITY) + Kratos.Array3([100.0, 100.0, 100.0]), 12)

        c += 200.0
        c.Evaluate(Kratos.ACCELERATION)
        for node in c.GetContainer():
            self.assertVectorAlmostEqual(self._GetValue(node, Kratos.ACCELERATION), self._GetValue(node, Kratos.VELOCITY) + Kratos.Array3([300.0, 300.0, 300.0]), 12)

        a = self._GetContainerExpression()
        b = self._GetContainerExpression()

        a.Read(Kratos.PRESSURE)
        b.Read(Kratos.PRESSURE)

        c = a + b
        c += b

        c.Evaluate(Kratos.DENSITY)
        for node in c.GetContainer():
            self.assertEqual(self._GetValue(node, Kratos.DENSITY), self._GetValue(node, Kratos.PRESSURE) * 3, 12)

        c = a + 100.0
        c.Evaluate(Kratos.DENSITY)
        for node in c.GetContainer():
            self.assertEqual(self._GetValue(node, Kratos.DENSITY), self._GetValue(node, Kratos.PRESSURE) + 100.0, 12)

        c += 100.0
        c.Evaluate(Kratos.DENSITY)
        for node in c.GetContainer():
            self.assertEqual(self._GetValue(node, Kratos.DENSITY), self._GetValue(node, Kratos.PRESSURE) + 200.0, 12)

    def test_ContainerExpressionMultiplyAndSubstract(self):
        a = self._GetContainerExpression()
        b = self._GetContainerExpression()

        a.Read(Kratos.VELOCITY)
        b.Read(Kratos.VELOCITY)

        c = a * 4 - b
        c *= 2
        c -= a

        c.Evaluate(Kratos.ACCELERATION)
        for node in c.GetContainer():
            self.assertVectorAlmostEqual(self._GetValue(node, Kratos.ACCELERATION), self._GetValue(node, Kratos.VELOCITY) * 5, 12)

        c = a - 100.0
        c.Evaluate(Kratos.ACCELERATION)
        for node in c.GetContainer():
            self.assertVectorAlmostEqual(self._GetValue(node, Kratos.ACCELERATION), self._GetValue(node, Kratos.VELOCITY) - Kratos.Array3([100.0, 100.0, 100.0]), 12)

        c -= 200.0
        c.Evaluate(Kratos.ACCELERATION)
        for node in c.GetContainer():
            self.assertVectorAlmostEqual(self._GetValue(node, Kratos.ACCELERATION), self._GetValue(node, Kratos.VELOCITY) - Kratos.Array3([300.0, 300.0, 300.0]), 12)

        a = self._GetContainerExpression()
        b = self._GetContainerExpression()

        a.Read(Kratos.PRESSURE)
        b.Read(Kratos.PRESSURE)

        c = a * 4 - b
        c *= 2
        c -= a

        c.Evaluate(Kratos.DENSITY)
        for node in c.GetContainer():
            self.assertEqual(self._GetValue(node, Kratos.DENSITY), self._GetValue(node, Kratos.PRESSURE) * 5, 12)

        c = a - 100.0
        c.Evaluate(Kratos.DENSITY)
        for node in c.GetContainer():
            self.assertEqual(self._GetValue(node, Kratos.DENSITY), self._GetValue(node, Kratos.PRESSURE) - 100.0, 12)

        c -= 100.0
        c.Evaluate(Kratos.DENSITY)
        for node in c.GetContainer():
            self.assertEqual(self._GetValue(node, Kratos.DENSITY), self._GetValue(node, Kratos.PRESSURE) - 200.0, 12)

        d = c * a
        d *= b
        d.Evaluate(Kratos.DENSITY)
        for node in c.GetContainer():
            self.assertEqual(self._GetValue(node, Kratos.DENSITY), (self._GetValue(node, Kratos.PRESSURE) - 200.0) * self._GetValue(node, Kratos.PRESSURE) ** 2, 12)

        a = self._GetContainerExpression()
        a.Read(Kratos.VELOCITY)

    def test_ContainerExpressionDivision(self):
        a = self._GetContainerExpression()
        a.Read(Kratos.VELOCITY)

        c = a / 2.0
        c /= 2.0

        c.Evaluate(Kratos.ACCELERATION)
        for node in c.GetContainer():
            self.assertVectorAlmostEqual(self._GetValue(node, Kratos.ACCELERATION), self._GetValue(node, Kratos.VELOCITY) / 4, 12)

        a = self._GetContainerExpression()
        a.Read(Kratos.PRESSURE)

        c = a / 2.0
        c /= 2.0

        c.Evaluate(Kratos.DENSITY)
        for node in c.GetContainer():
            self.assertEqual(self._GetValue(node, Kratos.DENSITY), self._GetValue(node, Kratos.PRESSURE) / 4, 12)

        d = c / a
        d /= (a * 2)
        d.Evaluate(Kratos.DENSITY)
        for node in c.GetContainer():
            self.assertEqual(self._GetValue(node, Kratos.DENSITY), 0.5 * ((self._GetValue(node, Kratos.PRESSURE) / 4) / self._GetValue(node, Kratos.PRESSURE)) / self._GetValue(node, Kratos.PRESSURE) , 12)

        a = self._GetContainerExpression()
        a.Read(Kratos.VELOCITY)

    def test_ContainerExpressionPow(self):
        a = self._GetContainerExpression()
        a.Read(Kratos.VELOCITY)

        c = a ** 2.0
        c **= 2.0

        c.Evaluate(Kratos.ACCELERATION)
        for node in c.GetContainer():
            ref_value = self._GetValue(node, Kratos.VELOCITY)
            self.assertVectorAlmostEqual(self._GetValue(node, Kratos.ACCELERATION), Kratos.Array3([ref_value[0]**4, ref_value[1]**4, ref_value[2]**4]), 12)

        a = self._GetContainerExpression()
        a.Read(Kratos.PRESSURE)

        c = a ** 2.0
        c **= 2.0

        c.Evaluate(Kratos.DENSITY)
        for node in c.GetContainer():
            self.assertEqual(self._GetValue(node, Kratos.DENSITY), self._GetValue(node, Kratos.PRESSURE) ** 4, 12)

    def test_ContainerExpressionNeg(self):
        a = self._GetContainerExpression()
        a.Read(Kratos.VELOCITY)

        c = -a

        c.Evaluate(Kratos.ACCELERATION)
        for node in c.GetContainer():
            self.assertVectorAlmostEqual(self._GetValue(node, Kratos.ACCELERATION), self._GetValue(node, Kratos.VELOCITY) * (-1.0), 12)

        a = self._GetContainerExpression()
        a.Read(Kratos.PRESSURE)

        c = -a

        c.Evaluate(Kratos.DENSITY)
        for node in c.GetContainer():
            self.assertEqual(self._GetValue(node, Kratos.DENSITY), self._GetValue(node, Kratos.PRESSURE) * (-1.0), 12)

    def test_SetDataForContainerVariable(self):
        a = self._GetContainerExpression()
        a.SetData(Kratos.Array3([1, 2, 3]))
        a.Evaluate(Kratos.ACCELERATION)
        for node in a.GetContainer():
            self.assertVectorAlmostEqual(self._GetValue(node, Kratos.ACCELERATION), Kratos.Array3([1, 2, 3]), 12)

        a = self._GetContainerExpression()
        a.SetData(10)
        a.Evaluate(Kratos.DENSITY)
        for node in a.GetContainer():
            self.assertEqual(self._GetValue(node, Kratos.DENSITY), 10)

    def test_ReadEvaluate(self) -> None:
        a = self._GetContainerExpression()
        a.Read(Kratos.VELOCITY)
        a.Evaluate(Kratos.ACCELERATION)
        for node in a.GetContainer():
            self.assertVectorAlmostEqual(self._GetValue(node, Kratos.ACCELERATION),
                                         self._GetValue(node, Kratos.VELOCITY),
                                         12)

    def test_Clone(self):
        a = self._GetContainerExpression()

        a.Read(Kratos.VELOCITY)

        b = a.Clone()
        b.SetData(Kratos.Array3([10, 11, 12]))

        a.Evaluate(Kratos.ACCELERATION)
        for node in a.GetContainer():
            self.assertVectorAlmostEqual(self._GetValue(node, Kratos.ACCELERATION), self._GetValue(node, Kratos.VELOCITY), 12)

        b.Evaluate(Kratos.ACCELERATION)
        for node in a.GetContainer():
            self.assertVectorAlmostEqual(self._GetValue(node, Kratos.ACCELERATION), Kratos.Array3([10, 11, 12]), 12)

        a = self._GetContainerExpression()

        a.Read(Kratos.PRESSURE)

        b = a.Clone()
        b.SetData(12)

        a.Evaluate(Kratos.DENSITY)
        for node in a.GetContainer():
            self.assertEqual(self._GetValue(node, Kratos.DENSITY), self._GetValue(node, Kratos.PRESSURE), 12)

        b.Evaluate(Kratos.DENSITY)
        for node in a.GetContainer():
            self.assertEqual(self._GetValue(node, Kratos.DENSITY), 12, 12)

    def test_GetContainer(self):
        a = self._GetContainerExpression()
        self.assertEqual(self._GetContainer(), a.GetContainer())

    @abstractmethod
    def _GetContainerExpression(self):
        pass

    @abstractmethod
    def _GetContainer(self):
        pass

    @abstractmethod
    def _GetValue(self, entity, variable):
        pass

class TestConditionPropertiesExpression(kratos_unittest.TestCase, TestContainerExpression):
    @classmethod
    def setUpClass(cls):
        cls.CreateEntities()

    def test_CopyData(self):
        a = self._GetContainerExpression()
        b = Kratos.Expression.ConditionNonHistoricalExpression(self.model_part)

        a.Read(Kratos.VELOCITY)
        b.CopyFrom(a)

        b.Evaluate(Kratos.ACCELERATION)
        for node in b.GetContainer():
            self.assertVectorAlmostEqual(node.GetValue(Kratos.ACCELERATION), self._GetValue(node, Kratos.VELOCITY), 12)

        a = self._GetContainerExpression()

        a.Read(Kratos.PRESSURE)
        b.CopyFrom(a)

        b.Evaluate(Kratos.DENSITY)
        for node in b.GetContainer():
            self.assertEqual(node.GetValue(Kratos.DENSITY), self._GetValue(node, Kratos.PRESSURE), 12)

        b = Kratos.Expression.ConditionNonHistoricalExpression(a)
        b += 1
        b.Evaluate(Kratos.DENSITY)
        for node in b.GetContainer():
            self.assertEqual(node.GetValue(Kratos.DENSITY), self._GetValue(node, Kratos.PRESSURE) + 1, 12)

    def _GetContainerExpression(self):
        return KratosOA.ContainerExpression.ConditionPropertiesExpression(self.model_part)

    def _GetContainer(self):
        return self.model_part.GetCommunicator().LocalMesh().Conditions

    def _GetValue(self, entity, variable):
        return entity.Properties[variable]

class TestElementPropertiesExpression(kratos_unittest.TestCase, TestContainerExpression):
    @classmethod
    def setUpClass(cls):
        cls.CreateEntities()

    def test_CopyData(self):
        a = self._GetContainerExpression()
        b = Kratos.Expression.ElementNonHistoricalExpression(self.model_part)

        a.Read(Kratos.VELOCITY)
        b.CopyFrom(a)

        b.Evaluate(Kratos.ACCELERATION)
        for node in b.GetContainer():
            self.assertVectorAlmostEqual(node.GetValue(Kratos.ACCELERATION), self._GetValue(node, Kratos.VELOCITY), 12)

        a = self._GetContainerExpression()

        a.Read(Kratos.PRESSURE)
        b.CopyFrom(a)

        b.Evaluate(Kratos.DENSITY)
        for node in b.GetContainer():
            self.assertEqual(node.GetValue(Kratos.DENSITY), self._GetValue(node, Kratos.PRESSURE), 12)

        b = Kratos.Expression.ElementNonHistoricalExpression(a)
        b += 1
        b.Evaluate(Kratos.DENSITY)
        for node in b.GetContainer():
            self.assertEqual(node.GetValue(Kratos.DENSITY), self._GetValue(node, Kratos.PRESSURE) + 1, 12)

    def _GetContainerExpression(self):
        return KratosOA.ContainerExpression.ElementPropertiesExpression(self.model_part)

    def _GetContainer(self):
        return self.model_part.GetCommunicator().LocalMesh().Elements

    def _GetValue(self, entity, variable):
        return entity.Properties[variable]

if __name__ == "__main__":
    Kratos.Tester.SetVerbosity(Kratos.Tester.Verbosity.PROGRESS)  # TESTS_OUTPUTS
    kratos_unittest.main()