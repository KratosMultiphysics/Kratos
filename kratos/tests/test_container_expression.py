
import numpy
from abc import ABC
from abc import abstractmethod
from typing import Union
import KratosMultiphysics as Kratos
from KratosMultiphysics.testing.utilities import ReadModelPart

# Import KratosUnittest
import KratosMultiphysics.KratosUnittest as kratos_unittest

class TestContainerExpression(ABC):
    @classmethod
    def CreateEntities(cls):
        cls.model =  Kratos.Model()
        cls.model_part = cls.model.CreateModelPart("test")
        cls.model_part.AddNodalSolutionStepVariable(Kratos.DENSITY)
        cls.model_part.AddNodalSolutionStepVariable(Kratos.PRESSURE)
        cls.model_part.AddNodalSolutionStepVariable(Kratos.ACCELERATION)
        cls.model_part.AddNodalSolutionStepVariable(Kratos.VELOCITY)
        cls.model_part.AddNodalSolutionStepVariable(Kratos.INITIAL_STRAIN)
        cls.model_part.AddNodalSolutionStepVariable(Kratos.GREEN_LAGRANGE_STRAIN_TENSOR)
        cls.model_part.AddNodalSolutionStepVariable(Kratos.PENALTY)
        cls.model_part.AddNodalSolutionStepVariable(Kratos.PK2_STRESS_TENSOR)
        cls.model_part.AddNodalSolutionStepVariable(Kratos.STEP)
        cls.model_part.AddNodalSolutionStepVariable(Kratos.STATIONARY)
        with kratos_unittest.WorkFolderScope(".", __file__, True):
            ReadModelPart("auxiliar_files_for_python_unittest/mdpa_files/two_dim_symmetrical_square", cls.model_part)

        for node in cls.model_part.Nodes:
            id = node.Id
            node.SetSolutionStepValue(Kratos.VELOCITY, Kratos.Array3([id+3, id+4, id+5]))
            node.SetSolutionStepValue(Kratos.PRESSURE, id+3)
            node.SetSolutionStepValue(Kratos.INITIAL_STRAIN, Kratos.Vector([id+3, id+4, id+5, id+6, id+7, id+8]))
            node.SetSolutionStepValue(Kratos.GREEN_LAGRANGE_STRAIN_TENSOR, Kratos.Matrix([[id+3, id+4], [id+5, id+6]]))
            node.SetValue(Kratos.PRESSURE, id+3)
            node.SetValue(Kratos.VELOCITY, Kratos.Array3([id+3, id+4, id+5]))
            node.SetValue(Kratos.INITIAL_STRAIN, Kratos.Vector([id+3, id+4, id+5, id+6, id+7, id+8]))
            node.SetValue(Kratos.GREEN_LAGRANGE_STRAIN_TENSOR, Kratos.Matrix([[id+3, id+4], [id+5, id+6]]))

        for condition in cls.model_part.Conditions:
            id = condition.Id
            condition.SetValue(Kratos.PRESSURE, id+4)
            condition.SetValue(Kratos.VELOCITY, Kratos.Array3([id+5, id+6, id+7]))
            condition.SetValue(Kratos.INITIAL_STRAIN, Kratos.Vector([id+3, id+4, id+5, id+6, id+7, id+8]))
            condition.SetValue(Kratos.GREEN_LAGRANGE_STRAIN_TENSOR, Kratos.Matrix([[id+3, id+4], [id+5, id+6]]))

        for element in cls.model_part.Elements:
            id = element.Id
            element.SetValue(Kratos.PRESSURE, id+5)
            element.SetValue(Kratos.VELOCITY, Kratos.Array3([id+6, id+7, id+8]))
            element.SetValue(Kratos.INITIAL_STRAIN, Kratos.Vector([id+3, id+4, id+5, id+6, id+7, id+8]))
            element.SetValue(Kratos.GREEN_LAGRANGE_STRAIN_TENSOR, Kratos.Matrix([[id+3, id+4], [id+5, id+6]]))

    def test_ContainerExpressionAdd(self):
        a = self._GetSpecializedContainerExpression()
        b = self._GetSpecializedContainerExpression()

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

        a = self._GetSpecializedContainerExpression()
        b = self._GetSpecializedContainerExpression()

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
        a = self._GetSpecializedContainerExpression()
        b = self._GetSpecializedContainerExpression()

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

        a = self._GetSpecializedContainerExpression()
        b = self._GetSpecializedContainerExpression()

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

        a = self._GetSpecializedContainerExpression()
        a.Read(Kratos.VELOCITY)

    def test_ContainerExpressionDivision(self):
        a = self._GetSpecializedContainerExpression()
        a.Read(Kratos.VELOCITY)

        c = a / 2.0
        c /= 2.0

        c.Evaluate(Kratos.ACCELERATION)
        for node in c.GetContainer():
            self.assertVectorAlmostEqual(self._GetValue(node, Kratos.ACCELERATION), self._GetValue(node, Kratos.VELOCITY) / 4, 12)

        a = self._GetSpecializedContainerExpression()
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

        a = self._GetSpecializedContainerExpression()
        a.Read(Kratos.VELOCITY)

    def test_ContainerExpressionPow(self):
        a = self._GetSpecializedContainerExpression()
        a.Read(Kratos.VELOCITY)

        c = a ** 2.0
        c **= 2.0

        c.Evaluate(Kratos.ACCELERATION)
        for node in c.GetContainer():
            ref_value = self._GetValue(node, Kratos.VELOCITY)
            self.assertVectorAlmostEqual(self._GetValue(node, Kratos.ACCELERATION), Kratos.Array3([ref_value[0]**4, ref_value[1]**4, ref_value[2]**4]), 12)

        a = self._GetSpecializedContainerExpression()
        a.Read(Kratos.PRESSURE)

        c = a ** 2.0
        c **= 2.0

        c.Evaluate(Kratos.DENSITY)
        for node in c.GetContainer():
            self.assertEqual(self._GetValue(node, Kratos.DENSITY), self._GetValue(node, Kratos.PRESSURE) ** 4, 12)

    def test_ContainerExpressionNeg(self):
        a = self._GetSpecializedContainerExpression()
        a.Read(Kratos.VELOCITY)

        c = -a

        c.Evaluate(Kratos.ACCELERATION)
        for node in c.GetContainer():
            self.assertVectorAlmostEqual(self._GetValue(node, Kratos.ACCELERATION), self._GetValue(node, Kratos.VELOCITY) * (-1.0), 12)

        a = self._GetSpecializedContainerExpression()
        a.Read(Kratos.PRESSURE)

        c = -a

        c.Evaluate(Kratos.DENSITY)
        for node in c.GetContainer():
            self.assertEqual(self._GetValue(node, Kratos.DENSITY), self._GetValue(node, Kratos.PRESSURE) * (-1.0), 12)

    def test_SetData(self):
        a = self._GetSpecializedContainerExpression()
        a.SetData(Kratos.Array3([1, 2, 3]))
        a.Evaluate(Kratos.ACCELERATION)
        for node in a.GetContainer():
            self.assertVectorAlmostEqual(self._GetValue(node, Kratos.ACCELERATION), Kratos.Array3([1, 2, 3]), 12)

        a = self._GetSpecializedContainerExpression()
        a.SetData(10)
        a.Evaluate(Kratos.DENSITY)
        for node in a.GetContainer():
            self.assertEqual(self._GetValue(node, Kratos.DENSITY), 10)

    def test_Clone(self):
        a = self._GetSpecializedContainerExpression()

        a.Read(Kratos.VELOCITY)

        b = a.Clone()
        b.SetData(Kratos.Array3([10, 11, 12]))

        a.Evaluate(Kratos.ACCELERATION)
        for node in a.GetContainer():
            self.assertVectorAlmostEqual(self._GetValue(node, Kratos.ACCELERATION), self._GetValue(node, Kratos.VELOCITY), 12)

        b.Evaluate(Kratos.ACCELERATION)
        for node in a.GetContainer():
            self.assertVectorAlmostEqual(self._GetValue(node, Kratos.ACCELERATION), Kratos.Array3([10, 11, 12]), 12)

        a = self._GetSpecializedContainerExpression()

        a.Read(Kratos.PRESSURE)

        b = a.Clone()
        b.SetData(12)

        a.Evaluate(Kratos.DENSITY)
        for node in a.GetContainer():
            self.assertEqual(self._GetValue(node, Kratos.DENSITY), self._GetValue(node, Kratos.PRESSURE), 12)

        b.Evaluate(Kratos.DENSITY)
        for node in a.GetContainer():
            self.assertEqual(self._GetValue(node, Kratos.DENSITY), 12, 12)

    def test_WeightedProduct(self):
        a = self._GetSpecializedContainerExpression()
        b = self._GetSpecializedContainerExpression()

        a.Read(Kratos.VELOCITY)
        b.Read(Kratos.PRESSURE)

        c = a * b
        c.Evaluate(Kratos.ACCELERATION)

        for entity in a.GetContainer():
            self.assertVectorAlmostEqual(self._GetValue(entity, Kratos.VELOCITY) * self._GetValue(entity, Kratos.PRESSURE), self._GetValue(entity, Kratos.ACCELERATION), 12)

    def test_Vector(self):
        a = self._GetSpecializedContainerExpression()
        a.Read(Kratos.INITIAL_STRAIN)

        a *= 2
        a.Evaluate(Kratos.PENALTY)

        for entity in a.GetContainer():
            self.assertVectorAlmostEqual(self._GetValue(entity, Kratos.INITIAL_STRAIN) * 2.0, self._GetValue(entity, Kratos.PENALTY), 12)

    def test_Matrix(self):
        a = self._GetSpecializedContainerExpression()
        a.Read(Kratos.GREEN_LAGRANGE_STRAIN_TENSOR)

        a *= 2
        a.Evaluate(Kratos.PK2_STRESS_TENSOR)

        for entity in a.GetContainer():
            self.assertMatrixAlmostEqual(self._GetValue(entity, Kratos.GREEN_LAGRANGE_STRAIN_TENSOR) * 2.0, self._GetValue(entity, Kratos.PK2_STRESS_TENSOR), 12)

    def test_Scope(self):
        def func(a):
            b = a * 2
            c = b * 3
            return c * 4

        a = self._GetSpecializedContainerExpression()

        a.Read(Kratos.VELOCITY)
        b = func(a)
        b.Evaluate(Kratos.ACCELERATION)

        for entity in a.GetContainer():
            self.assertVectorAlmostEqual(self._GetValue(entity, Kratos.VELOCITY) * 24.0, self._GetValue(entity, Kratos.ACCELERATION), 12)

    def test_NumpyEvaluate(self):
        a = self._GetSpecializedContainerExpression()

        a.Read(Kratos.PRESSURE)
        numpy_array = a.Evaluate()
        for i, entity in enumerate(a.GetContainer()):
            self.assertEqual(self._GetValue(entity, Kratos.PRESSURE), numpy_array[i], 12)

        a.Read(Kratos.VELOCITY)
        numpy_array = a.Evaluate()
        for i, entity in enumerate(a.GetContainer()):
            self.assertVectorAlmostEqual(self._GetValue(entity, Kratos.VELOCITY), numpy_array[i, :], 12)

        a.Read(Kratos.INITIAL_STRAIN)
        numpy_array = a.Evaluate()
        for i, entity in enumerate(a.GetContainer()):
            self.assertVectorAlmostEqual(self._GetValue(entity, Kratos.INITIAL_STRAIN), numpy_array[i, :], 12)

        a.Read(Kratos.GREEN_LAGRANGE_STRAIN_TENSOR)
        numpy_array = a.Evaluate()
        for i_entity, entity in enumerate(a.GetContainer()):
            v = self._GetValue(entity, Kratos.GREEN_LAGRANGE_STRAIN_TENSOR)
            for i in range(v.Size1()):
                for j in range(v.Size2()):
                    self.assertEqual(v[i, j], numpy_array[i_entity, i, j])

    def test_NumpyRead(self):
        a = self._GetSpecializedContainerExpression()

        numpy_array = numpy.arange(0.0, len(a.GetContainer()))
        a.Read(numpy_array)
        a.Evaluate(Kratos.PRESSURE)

        for i, entity in enumerate(a.GetContainer()):
            v = self._GetValue(entity, Kratos.PRESSURE)
            self.assertEqual(v, numpy_array[i], 12)

        numpy_array = numpy.arange(0.0, len(a.GetContainer()) * 8).reshape(len(a.GetContainer()), 8)
        a.Read(numpy_array)
        a *= 2
        a.Evaluate(Kratos.PENALTY)
        for i, entity in enumerate(a.GetContainer()):
            self.assertVectorAlmostEqual(self._GetValue(entity, Kratos.PENALTY), numpy_array[i, :] * 2, 12)

        numpy_array = numpy.arange(0.0, len(a.GetContainer()) * 2 * 4).reshape(len(a.GetContainer()), 2, 4)
        a.Read(numpy_array)
        a *= 2
        a.Evaluate(Kratos.PK2_STRESS_TENSOR)
        for i_entity, entity in enumerate(a.GetContainer()):
            v = self._GetValue(entity, Kratos.PK2_STRESS_TENSOR)
            for i in range(v.Size1()):
                for j in range(v.Size2()):
                    self.assertEqual(v[i, j], numpy_array[i_entity, i, j] * 2)

    def test_NumpyMoveFrom(self):
        # double move check
        a = self._GetSpecializedContainerExpression()

        numpy_array = numpy.arange(0.0, len(a.GetContainer()))
        a.MoveFrom(numpy_array)
        a.Evaluate(Kratos.PRESSURE)

        for i, entity in enumerate(a.GetContainer()):
            v = self._GetValue(entity, Kratos.PRESSURE)
            self.assertEqual(v, numpy_array[i], 12)

        numpy_array[:] = numpy_array + 1

        a.Evaluate(Kratos.DENSITY)
        for i, entity in enumerate(a.GetContainer()):
            v = self._GetValue(entity, Kratos.DENSITY)
            self.assertEqual(v, numpy_array[i], 12)

        # int move check
        a = self._GetSpecializedContainerExpression()

        numpy_array = numpy.arange(0, len(a.GetContainer()), dtype=numpy.int32)
        a.MoveFrom(numpy_array)
        a.Evaluate(Kratos.STEP)

        for i, entity in enumerate(a.GetContainer()):
            v = self._GetValue(entity, Kratos.STEP)
            self.assertEqual(v, numpy_array[i], 12)

        numpy_array[:] = numpy_array + 1

        a.Evaluate(Kratos.STATIONARY)
        for i, entity in enumerate(a.GetContainer()):
            v = self._GetValue(entity, Kratos.STATIONARY)
            self.assertEqual(v, numpy_array[i], 12)

    def test_NumpyForbiddenCasts(self):
        a = self._GetSpecializedContainerExpression()

        with self.assertRaises(TypeError):
            numpy_array = numpy.arange(0, len(a.GetContainer()), dtype=numpy.float32)
            a.Read(numpy_array)

        with self.assertRaises(TypeError):
            numpy_array = numpy.arange(0, len(a.GetContainer()), dtype=numpy.int64)
            a.Read(numpy_array)

        with self.assertRaises(TypeError):
            numpy_array = numpy.arange(0, len(a.GetContainer()), dtype=numpy.float32)
            a.MoveFrom(numpy_array)

        with self.assertRaises(TypeError):
            numpy_array = numpy.arange(0, len(a.GetContainer()), dtype=numpy.int64)
            a.MoveFrom(numpy_array)


    def test_GetContainer(self):
        a = self._GetSpecializedContainerExpression()
        self.assertEqual(self._GetContainer(), a.GetContainer())

    @abstractmethod
    def _GetSpecializedContainerExpression(self) -> Union[Kratos.ContainerExpression.HistoricalExpression, Kratos.ContainerExpression.NodalNonHistoricalExpression, Kratos.ContainerExpression.ConditionNonHistoricalExpression, Kratos.ContainerExpression.ElementNonHistoricalExpression]:
        pass

    @abstractmethod
    def _GetContainer(self) -> Union[Kratos.NodesArray, Kratos.ConditionsArray, Kratos.ElementsArray]:
        pass

    @abstractmethod
    def _GetValue(self, entity, variable):
        pass

class TestHistoricalContainerExpression(kratos_unittest.TestCase, TestContainerExpression):
    @classmethod
    def setUpClass(cls):
        cls.CreateEntities()

    def test_CopyFrom(self):
        a = self._GetSpecializedContainerExpression()
        b = Kratos.ContainerExpression.NodalNonHistoricalExpression(self.model_part)

        a.Read(Kratos.VELOCITY)
        b.CopyFrom(a)

        b.Evaluate(Kratos.ACCELERATION)
        for node in b.GetContainer():
            self.assertVectorAlmostEqual(node.GetValue(Kratos.ACCELERATION), self._GetValue(node, Kratos.VELOCITY), 12)

        a = self._GetSpecializedContainerExpression()

        a.Read(Kratos.PRESSURE)
        b.CopyFrom(a)

        b.Evaluate(Kratos.DENSITY)
        for node in b.GetContainer():
            self.assertEqual(node.GetValue(Kratos.DENSITY), self._GetValue(node, Kratos.PRESSURE), 12)

        b = Kratos.ContainerExpression.NodalNonHistoricalExpression(a)
        b += 1
        b.Evaluate(Kratos.DENSITY)
        for node in b.GetContainer():
            self.assertEqual(node.GetValue(Kratos.DENSITY), self._GetValue(node, Kratos.PRESSURE) + 1, 12)

    def _GetSpecializedContainerExpression(self):
        return Kratos.ContainerExpression.HistoricalExpression(self.model_part)

    def _GetContainer(self):
        return self.model_part.GetCommunicator().LocalMesh().Nodes

    def _GetValue(self, entity, variable):
        return entity.GetSolutionStepValue(variable)

class TestNodalContainerExpression(kratos_unittest.TestCase, TestContainerExpression):
    @classmethod
    def setUpClass(cls):
        cls.CreateEntities()

    def test_CopyFrom(self):
        a = self._GetSpecializedContainerExpression()
        b = Kratos.ContainerExpression.HistoricalExpression(self.model_part)

        a.Read(Kratos.VELOCITY)
        b.CopyFrom(a)

        b.Evaluate(Kratos.ACCELERATION)
        for node in b.GetContainer():
            self.assertVectorAlmostEqual(node.GetSolutionStepValue(Kratos.ACCELERATION), self._GetValue(node, Kratos.VELOCITY), 12)

        a = self._GetSpecializedContainerExpression()

        a.Read(Kratos.PRESSURE)
        b.CopyFrom(a)

        b.Evaluate(Kratos.DENSITY)
        for node in b.GetContainer():
            self.assertEqual(node.GetSolutionStepValue(Kratos.DENSITY), self._GetValue(node, Kratos.PRESSURE), 12)

        b = Kratos.ContainerExpression.HistoricalExpression(a)
        b += 1
        b.Evaluate(Kratos.DENSITY)
        for node in b.GetContainer():
            self.assertEqual(node.GetSolutionStepValue(Kratos.DENSITY), self._GetValue(node, Kratos.PRESSURE) + 1, 12)

    def _GetSpecializedContainerExpression(self):
        return Kratos.ContainerExpression.NodalNonHistoricalExpression(self.model_part)

    def _GetContainer(self):
        return self.model_part.GetCommunicator().LocalMesh().Nodes

    def _GetValue(self, entity, variable):
        return entity.GetValue(variable)

class TestConditionContainerExpression(kratos_unittest.TestCase, TestContainerExpression):
    @classmethod
    def setUpClass(cls):
        cls.CreateEntities()

    def _GetSpecializedContainerExpression(self):
        return Kratos.ContainerExpression.ConditionNonHistoricalExpression(self.model_part)

    def _GetContainer(self):
        return self.model_part.GetCommunicator().LocalMesh().Conditions

    def _GetValue(self, entity, variable):
        return entity.GetValue(variable)

class TestElementContainerExpression(kratos_unittest.TestCase, TestContainerExpression):
    @classmethod
    def setUpClass(cls):
        cls.CreateEntities()

    def _GetSpecializedContainerExpression(self):
        return Kratos.ContainerExpression.ElementNonHistoricalExpression(self.model_part)

    def _GetContainer(self):
        return self.model_part.GetCommunicator().LocalMesh().Elements

    def _GetValue(self, entity, variable):
        return entity.GetValue(variable)

if __name__ == "__main__":
    Kratos.Tester.SetVerbosity(Kratos.Tester.Verbosity.PROGRESS)  # TESTS_OUTPUTS
    kratos_unittest.main()
