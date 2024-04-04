
import numpy
from abc import ABC
from abc import abstractmethod
from typing import Union
import KratosMultiphysics as Kratos
from KratosMultiphysics.testing.utilities import ReadModelPart

# Import KratosUnittest
import KratosMultiphysics.KratosUnittest as kratos_unittest

class TestContainerExpression(ABC):
    ExpressionUnionType = Union[Kratos.Expression.NodalExpression, Kratos.Expression.ConditionExpression, Kratos.Expression.ElementExpression]
    @classmethod
    def CreateEntities(cls):
        cls.model =  Kratos.Model()
        cls.model_part = cls.model.CreateModelPart("test")
        cls.model_part.AddNodalSolutionStepVariable(Kratos.DENSITY)
        cls.model_part.AddNodalSolutionStepVariable(Kratos.PRESSURE)
        cls.model_part.AddNodalSolutionStepVariable(Kratos.TEMPERATURE)
        cls.model_part.AddNodalSolutionStepVariable(Kratos.ACCELERATION)
        cls.model_part.AddNodalSolutionStepVariable(Kratos.VELOCITY)
        cls.model_part.AddNodalSolutionStepVariable(Kratos.REACTION)
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

        cls.data_comm: Kratos.DataCommunicator = cls.model_part.GetCommunicator().GetDataCommunicator()

    def test_ContainerExpressionAdd(self):
        a = self._GetContainerExpression()
        b = self._GetContainerExpression()

        self._Read(a, Kratos.VELOCITY)
        self._Read(b, Kratos.VELOCITY)

        c = a + b
        c += b

        self._Evaluate(c, Kratos.ACCELERATION)
        self._Evaluate(a, Kratos.REACTION)
        for node in c.GetContainer():
            self.assertVectorAlmostEqual(self._GetValue(node, Kratos.ACCELERATION), self._GetValue(node, Kratos.VELOCITY) * 3, 12)
            self.assertVectorAlmostEqual(self._GetValue(node, Kratos.REACTION), self._GetValue(node, Kratos.VELOCITY), 12)

        c = a + 100.0
        self._Evaluate(c, Kratos.ACCELERATION)
        for node in c.GetContainer():
            self.assertVectorAlmostEqual(self._GetValue(node, Kratos.ACCELERATION), self._GetValue(node, Kratos.VELOCITY) + Kratos.Array3([100.0, 100.0, 100.0]), 12)

        c += 200.0
        self._Evaluate(c, Kratos.ACCELERATION)
        for node in c.GetContainer():
            self.assertVectorAlmostEqual(self._GetValue(node, Kratos.ACCELERATION), self._GetValue(node, Kratos.VELOCITY) + Kratos.Array3([300.0, 300.0, 300.0]), 12)

        a = self._GetContainerExpression()
        b = self._GetContainerExpression()

        self._Read(a, Kratos.PRESSURE)
        self._Read(b, Kratos.PRESSURE)

        c = a + b
        c += b

        self._Evaluate(c, Kratos.DENSITY)
        for node in c.GetContainer():
            self.assertEqual(self._GetValue(node, Kratos.DENSITY), self._GetValue(node, Kratos.PRESSURE) * 3, 12)

        c = a + 100.0
        self._Evaluate(c, Kratos.DENSITY)
        for node in c.GetContainer():
            self.assertEqual(self._GetValue(node, Kratos.DENSITY), self._GetValue(node, Kratos.PRESSURE) + 100.0, 12)

        c += 100.0
        self._Evaluate(c, Kratos.DENSITY)
        for node in c.GetContainer():
            self.assertEqual(self._GetValue(node, Kratos.DENSITY), self._GetValue(node, Kratos.PRESSURE) + 200.0, 12)

    def test_ContainerExpressionMultiplyAndSubstract(self):
        a = self._GetContainerExpression()
        b = self._GetContainerExpression()

        self._Read(a, Kratos.VELOCITY)
        self._Read(b, Kratos.VELOCITY)

        c = a * 4 - b
        c *= 2
        c -= a

        self._Evaluate(c, Kratos.ACCELERATION)
        self._Evaluate(a, Kratos.REACTION)
        for node in c.GetContainer():
            self.assertVectorAlmostEqual(self._GetValue(node, Kratos.ACCELERATION), self._GetValue(node, Kratos.VELOCITY) * 5, 12)
            self.assertVectorAlmostEqual(self._GetValue(node, Kratos.REACTION), self._GetValue(node, Kratos.VELOCITY), 12)

        c = a - 100.0
        self._Evaluate(c, Kratos.ACCELERATION)
        for node in c.GetContainer():
            self.assertVectorAlmostEqual(self._GetValue(node, Kratos.ACCELERATION), self._GetValue(node, Kratos.VELOCITY) - Kratos.Array3([100.0, 100.0, 100.0]), 12)

        c -= 200.0
        self._Evaluate(c, Kratos.ACCELERATION)
        for node in c.GetContainer():
            self.assertVectorAlmostEqual(self._GetValue(node, Kratos.ACCELERATION), self._GetValue(node, Kratos.VELOCITY) - Kratos.Array3([300.0, 300.0, 300.0]), 12)

        a = self._GetContainerExpression()
        b = self._GetContainerExpression()

        self._Read(a, Kratos.PRESSURE)
        self._Read(b, Kratos.PRESSURE)

        c = a * 4 - b
        c *= 2
        c -= a

        self._Evaluate(c, Kratos.DENSITY)
        for node in c.GetContainer():
            self.assertEqual(self._GetValue(node, Kratos.DENSITY), self._GetValue(node, Kratos.PRESSURE) * 5, 12)

        c = a - 100.0
        self._Evaluate(c, Kratos.DENSITY)
        for node in c.GetContainer():
            self.assertEqual(self._GetValue(node, Kratos.DENSITY), self._GetValue(node, Kratos.PRESSURE) - 100.0, 12)

        c -= 100.0
        self._Evaluate(c, Kratos.DENSITY)
        for node in c.GetContainer():
            self.assertEqual(self._GetValue(node, Kratos.DENSITY), self._GetValue(node, Kratos.PRESSURE) - 200.0, 12)

        d = c * a
        d *= b
        self._Evaluate(d, Kratos.DENSITY)
        for node in c.GetContainer():
            self.assertEqual(self._GetValue(node, Kratos.DENSITY), (self._GetValue(node, Kratos.PRESSURE) - 200.0) * self._GetValue(node, Kratos.PRESSURE) ** 2, 12)

        a = self._GetContainerExpression()
        self._Read(a, Kratos.VELOCITY)

    def test_ContainerExpressionDivision(self):
        a = self._GetContainerExpression()
        self._Read(a, Kratos.VELOCITY)

        c = a / 2.0
        c /= 2.0

        self._Evaluate(c, Kratos.ACCELERATION)
        self._Evaluate(a, Kratos.REACTION)
        for node in c.GetContainer():
            self.assertVectorAlmostEqual(self._GetValue(node, Kratos.ACCELERATION), self._GetValue(node, Kratos.VELOCITY) / 4, 12)
            self.assertVectorAlmostEqual(self._GetValue(node, Kratos.REACTION), self._GetValue(node, Kratos.VELOCITY), 12)

        a = self._GetContainerExpression()
        self._Read(a, Kratos.PRESSURE)

        c = a / 2.0
        c /= 2.0

        self._Evaluate(c, Kratos.DENSITY)
        self._Evaluate(a, Kratos.TEMPERATURE)
        for node in c.GetContainer():
            self.assertEqual(self._GetValue(node, Kratos.DENSITY), self._GetValue(node, Kratos.PRESSURE) / 4, 12)
            self.assertEqual(self._GetValue(node, Kratos.TEMPERATURE), self._GetValue(node, Kratos.PRESSURE), 12)

        d = c / a
        d /= (a * 2)
        self._Evaluate(d, Kratos.DENSITY)
        for node in c.GetContainer():
            self.assertEqual(self._GetValue(node, Kratos.DENSITY), 0.5 * ((self._GetValue(node, Kratos.PRESSURE) / 4) / self._GetValue(node, Kratos.PRESSURE)) / self._GetValue(node, Kratos.PRESSURE) , 12)

    def test_ContainerExpressionPow(self):
        a = self._GetContainerExpression()
        self._Read(a, Kratos.VELOCITY)

        c = a ** 2.0
        c **= 2.0

        self._Evaluate(c, Kratos.ACCELERATION)
        self._Evaluate(a, Kratos.REACTION)
        for node in c.GetContainer():
            ref_value = self._GetValue(node, Kratos.VELOCITY)
            self.assertVectorAlmostEqual(self._GetValue(node, Kratos.ACCELERATION), Kratos.Array3([ref_value[0]**4, ref_value[1]**4, ref_value[2]**4]), 12)
            self.assertVectorAlmostEqual(self._GetValue(node, Kratos.REACTION), ref_value, 12)

        a = self._GetContainerExpression()
        self._Read(a, Kratos.PRESSURE)

        c = a ** 2.0
        c **= 2.0

        self._Evaluate(c, Kratos.DENSITY)
        self._Evaluate(a, Kratos.TEMPERATURE)
        for node in c.GetContainer():
            self.assertEqual(self._GetValue(node, Kratos.DENSITY), self._GetValue(node, Kratos.PRESSURE) ** 4, 12)
            self.assertEqual(self._GetValue(node, Kratos.TEMPERATURE), self._GetValue(node, Kratos.PRESSURE), 12)

    def test_ContainerExpressionNeg(self):
        a = self._GetContainerExpression()
        self._Read(a, Kratos.VELOCITY)

        c = -a

        self._Evaluate(c, Kratos.ACCELERATION)
        self._Evaluate(a, Kratos.REACTION)
        for node in c.GetContainer():
            self.assertVectorAlmostEqual(self._GetValue(node, Kratos.ACCELERATION), self._GetValue(node, Kratos.VELOCITY) * (-1.0), 12)
            self.assertVectorAlmostEqual(self._GetValue(node, Kratos.REACTION), self._GetValue(node, Kratos.VELOCITY), 12)

        a = self._GetContainerExpression()
        self._Read(a, Kratos.PRESSURE)

        c = -a

        self._Evaluate(c, Kratos.DENSITY)
        self._Evaluate(a, Kratos.TEMPERATURE)
        for node in c.GetContainer():
            self.assertEqual(self._GetValue(node, Kratos.DENSITY), self._GetValue(node, Kratos.PRESSURE) * (-1.0), 12)
            self.assertEqual(self._GetValue(node, Kratos.TEMPERATURE), self._GetValue(node, Kratos.PRESSURE), 12)

    def test_SetData(self):
        a = self._GetContainerExpression()
        Kratos.Expression.LiteralExpressionIO.SetData(a, Kratos.Array3([1, 2, 3]))
        self._Evaluate(a, Kratos.ACCELERATION)
        for node in a.GetContainer():
            self.assertVectorAlmostEqual(self._GetValue(node, Kratos.ACCELERATION), Kratos.Array3([1, 2, 3]), 12)

        a = self._GetContainerExpression()
        Kratos.Expression.LiteralExpressionIO.SetData(a, 10)
        self._Evaluate(a, Kratos.DENSITY)
        for node in a.GetContainer():
            self.assertEqual(self._GetValue(node, Kratos.DENSITY), 10)

    def test_Clone(self):
        a = self._GetContainerExpression()

        self._Read(a, Kratos.VELOCITY)

        b = a.Clone()
        Kratos.Expression.LiteralExpressionIO.SetData(b, Kratos.Array3([10, 11, 12]))

        self._Evaluate(a, Kratos.ACCELERATION)
        for node in a.GetContainer():
            self.assertVectorAlmostEqual(self._GetValue(node, Kratos.ACCELERATION), self._GetValue(node, Kratos.VELOCITY), 12)

        self._Evaluate(b, Kratos.ACCELERATION)
        for node in a.GetContainer():
            self.assertVectorAlmostEqual(self._GetValue(node, Kratos.ACCELERATION), Kratos.Array3([10, 11, 12]), 12)

        a = self._GetContainerExpression()

        self._Read(a, Kratos.PRESSURE)

        b = a.Clone()
        Kratos.Expression.LiteralExpressionIO.SetData(b, 12)

        self._Evaluate(a, Kratos.DENSITY)
        for node in a.GetContainer():
            self.assertEqual(self._GetValue(node, Kratos.DENSITY), self._GetValue(node, Kratos.PRESSURE), 12)

        self._Evaluate(b, Kratos.DENSITY)
        for node in a.GetContainer():
            self.assertEqual(self._GetValue(node, Kratos.DENSITY), 12, 12)

    def test_WeightedProduct(self):
        a = self._GetContainerExpression()
        b = self._GetContainerExpression()

        self._Read(a, Kratos.VELOCITY)
        self._Read(b, Kratos.PRESSURE)

        c = Kratos.Expression.Utils.Scale(a, b)
        self._Evaluate(c, Kratos.ACCELERATION)

        for entity in a.GetContainer():
            self.assertVectorAlmostEqual(self._GetValue(entity, Kratos.VELOCITY) * self._GetValue(entity, Kratos.PRESSURE), self._GetValue(entity, Kratos.ACCELERATION), 12)

    def test_Vector(self):
        a = self._GetContainerExpression()
        self._Read(a, Kratos.INITIAL_STRAIN)

        a *= 2
        self._Evaluate(a, Kratos.PENALTY)

        for entity in a.GetContainer():
            self.assertVectorAlmostEqual(self._GetValue(entity, Kratos.INITIAL_STRAIN) * 2.0, self._GetValue(entity, Kratos.PENALTY), 12)

    def test_Matrix(self):
        a = self._GetContainerExpression()
        self._Read(a, Kratos.GREEN_LAGRANGE_STRAIN_TENSOR)

        a *= 2
        self._Evaluate(a, Kratos.PK2_STRESS_TENSOR)

        for entity in a.GetContainer():
            self.assertMatrixAlmostEqual(self._GetValue(entity, Kratos.GREEN_LAGRANGE_STRAIN_TENSOR) * 2.0, self._GetValue(entity, Kratos.PK2_STRESS_TENSOR), 12)

    def test_Scope(self):
        def func(a):
            b = a * 2
            c = b * 3
            return c * 4

        a = self._GetContainerExpression()

        self._Read(a, Kratos.VELOCITY)
        b = func(a)
        self._Evaluate(b, Kratos.ACCELERATION)

        for entity in a.GetContainer():
            self.assertVectorAlmostEqual(self._GetValue(entity, Kratos.VELOCITY) * 24.0, self._GetValue(entity, Kratos.ACCELERATION), 12)

    def test_NumpyEvaluate(self):
        """@todo: Rmove (replaced by test_NumpyEvaluateContainerExpression)"""
        a = self._GetContainerExpression()

        self._Read(a, Kratos.PRESSURE)
        numpy_array = a.Evaluate()
        for i, entity in enumerate(a.GetContainer()):
            self.assertEqual(self._GetValue(entity, Kratos.PRESSURE), numpy_array[i], 12)

        self._Read(a, Kratos.VELOCITY)
        numpy_array = a.Evaluate()
        for i, entity in enumerate(a.GetContainer()):
            self.assertVectorAlmostEqual(self._GetValue(entity, Kratos.VELOCITY), numpy_array[i, :], 12)

        self._Read(a, Kratos.INITIAL_STRAIN)
        numpy_array = a.Evaluate()
        for i, entity in enumerate(a.GetContainer()):
            self.assertVectorAlmostEqual(self._GetValue(entity, Kratos.INITIAL_STRAIN), numpy_array[i, :], 12)

        self._Read(a, Kratos.GREEN_LAGRANGE_STRAIN_TENSOR)
        numpy_array = a.Evaluate()
        for i_entity, entity in enumerate(a.GetContainer()):
            v = self._GetValue(entity, Kratos.GREEN_LAGRANGE_STRAIN_TENSOR)
            for i in range(v.Size1()):
                for j in range(v.Size2()):
                    self.assertEqual(v[i, j], numpy_array[i_entity, i, j])

    def test_NumpyEvaluateContainerExpression(self):
        a = self._GetContainerExpression()

        self._Read(a, Kratos.PRESSURE)
        numpy_array = a.Evaluate()
        for i, entity in enumerate(a.GetContainer()):
            self.assertEqual(self._GetValue(entity, Kratos.PRESSURE), numpy_array[i], 12)

        self._Read(a, Kratos.VELOCITY)
        numpy_array = a.Evaluate()

        for i, entity in enumerate(a.GetContainer()):
            self.assertVectorAlmostEqual(self._GetValue(entity, Kratos.VELOCITY), numpy_array[i, :], 12)

        self._Read(a, Kratos.INITIAL_STRAIN)
        numpy_array = a.Evaluate()

        for i, entity in enumerate(a.GetContainer()):
            self.assertVectorAlmostEqual(self._GetValue(entity, Kratos.INITIAL_STRAIN), numpy_array[i, :], 12)

        self._Read(a, Kratos.GREEN_LAGRANGE_STRAIN_TENSOR)
        numpy_array = a.Evaluate()

        for i_entity, entity in enumerate(a.GetContainer()):
            v = self._GetValue(entity, Kratos.GREEN_LAGRANGE_STRAIN_TENSOR)
            for i in range(v.Size1()):
                for j in range(v.Size2()):
                    self.assertEqual(v[i, j], numpy_array[i_entity, i, j])

    def test_NumpyRead(self):
        """@todo Remove (replaced by test_NumpyReadContainerExpression)"""
        a = self._GetContainerExpression()

        numpy_array = numpy.arange(0.0, len(a.GetContainer()))
        Kratos.Expression.CArrayExpressionIO.Read(a, numpy_array)
        self._Evaluate(a, Kratos.PRESSURE)

        for i, entity in enumerate(a.GetContainer()):
            v = self._GetValue(entity, Kratos.PRESSURE)
            self.assertEqual(v, numpy_array[i], 12)

        numpy_array = numpy.arange(0.0, len(a.GetContainer()) * 8).reshape(len(a.GetContainer()), 8)
        Kratos.Expression.CArrayExpressionIO.Read(a, numpy_array)
        a *= 2
        self._Evaluate(a, Kratos.PENALTY)
        for i, entity in enumerate(a.GetContainer()):
            self.assertVectorAlmostEqual(self._GetValue(entity, Kratos.PENALTY), numpy_array[i, :] * 2, 12)

        numpy_array = numpy.arange(0.0, len(a.GetContainer()) * 2 * 4).reshape(len(a.GetContainer()), 2, 4)
        Kratos.Expression.CArrayExpressionIO.Read(a, numpy_array)
        a *= 2
        self._Evaluate(a, Kratos.PK2_STRESS_TENSOR)
        for i_entity, entity in enumerate(a.GetContainer()):
            v = self._GetValue(entity, Kratos.PK2_STRESS_TENSOR)
            for i in range(v.Size1()):
                for j in range(v.Size2()):
                    self.assertEqual(v[i, j], numpy_array[i_entity, i, j] * 2)

    def test_NumpyReadContainerExpression(self):
        a = self._GetContainerExpression()

        numpy_array = numpy.arange(0.0, len(a.GetContainer()))
        Kratos.Expression.CArrayExpressionIO.Read(a, numpy_array)
        self._Evaluate(a, Kratos.PRESSURE)
        for i, entity in enumerate(a.GetContainer()):
            v = self._GetValue(entity, Kratos.PRESSURE)
            self.assertEqual(v, numpy_array[i], 12)

        numpy_array = numpy.arange(0.0, len(a.GetContainer()) * 8).reshape(len(a.GetContainer()), 8)
        Kratos.Expression.CArrayExpressionIO.Read(a, numpy_array)
        a *= 2
        self._Evaluate(a, Kratos.PENALTY)
        for i, entity in enumerate(a.GetContainer()):
            self.assertVectorAlmostEqual(self._GetValue(entity, Kratos.PENALTY), numpy_array[i, :] * 2, 12)

        numpy_array = numpy.arange(0.0, len(a.GetContainer()) * 2 * 4).reshape(len(a.GetContainer()), 2, 4)
        Kratos.Expression.CArrayExpressionIO.Read(a, numpy_array)
        a *= 2
        self._Evaluate(a, Kratos.PK2_STRESS_TENSOR)
        for i_entity, entity in enumerate(a.GetContainer()):
            v = self._GetValue(entity, Kratos.PK2_STRESS_TENSOR)
            for i in range(v.Size1()):
                for j in range(v.Size2()):
                    self.assertEqual(v[i, j], numpy_array[i_entity, i, j] * 2)

    def test_NumpyMoveFrom(self):
        # double move check
        a = self._GetContainerExpression()

        numpy_array = numpy.arange(0.0, len(a.GetContainer()))
        Kratos.Expression.CArrayExpressionIO.Move(a, numpy_array)
        self._Evaluate(a, Kratos.PRESSURE)

        for i, entity in enumerate(a.GetContainer()):
            v = self._GetValue(entity, Kratos.PRESSURE)
            self.assertEqual(v, numpy_array[i], 12)

        numpy_array[:] = numpy_array + 1

        self._Evaluate(a, Kratos.DENSITY)
        for i, entity in enumerate(a.GetContainer()):
            v = self._GetValue(entity, Kratos.DENSITY)
            self.assertEqual(v, numpy_array[i], 12)

        # int move check
        a = self._GetContainerExpression()

        numpy_array = numpy.arange(0, len(a.GetContainer()), dtype=numpy.int32)
        Kratos.Expression.CArrayExpressionIO.Move(a, numpy_array)
        self._Evaluate(a, Kratos.STEP)

        for i, entity in enumerate(a.GetContainer()):
            v = self._GetValue(entity, Kratos.STEP)
            self.assertEqual(v, numpy_array[i], 12)

        numpy_array[:] = numpy_array + 1

        self._Evaluate(a, Kratos.STATIONARY)
        for i, entity in enumerate(a.GetContainer()):
            v = self._GetValue(entity, Kratos.STATIONARY)
            self.assertEqual(v, numpy_array[i], 12)

    def test_NumpyForbiddenCasts(self):
        a = self._GetContainerExpression()

        with self.assertRaises(TypeError):
            numpy_array = numpy.arange(0, len(a.GetContainer()), dtype=numpy.float32)
            self._Read(a, numpy_array)

        with self.assertRaises(TypeError):
            numpy_array = numpy.arange(0, len(a.GetContainer()), dtype=numpy.int64)
            self._Read(a, numpy_array)

        with self.assertRaises(TypeError):
            numpy_array = numpy.arange(0, len(a.GetContainer()), dtype=numpy.float32)
            Kratos.Expression.CArrayExpressionIO.Move(a, numpy_array)

        with self.assertRaises(TypeError):
            numpy_array = numpy.arange(0, len(a.GetContainer()), dtype=numpy.int64)
            Kratos.Expression.CArrayExpressionIO.Move(a, numpy_array)

    def test_Slice(self):
        a = self._GetContainerExpression()
        self._Read(a, Kratos.INITIAL_STRAIN)

        sliced = Kratos.Expression.Utils.Slice(a, 2, 3)
        sliced *= 2.0
        self._Evaluate(sliced, Kratos.ACCELERATION)

        for entity in a.GetContainer():
            original_value = self._GetValue(entity, Kratos.INITIAL_STRAIN)
            new_array = Kratos.Array3([original_value[2], original_value[3], original_value[4]])
            self.assertVectorAlmostEqual(self._GetValue(entity, Kratos.ACCELERATION), new_array * 2, 12)

    def test_Reshape(self):
        a = self._GetContainerExpression()
        self._Read(a, Kratos.INITIAL_STRAIN)

        reshaped = Kratos.Expression.Utils.Reshape(a, [2, 3])
        reshaped *= 2.0
        self._Evaluate(reshaped, Kratos.PK2_STRESS_TENSOR)

        for entity in a.GetContainer():
            original_value = self._GetValue(entity, Kratos.INITIAL_STRAIN)
            new_matrix = Kratos.Matrix(2, 3)
            new_matrix[0, 0] = original_value[0]
            new_matrix[0, 1] = original_value[1]
            new_matrix[0, 2] = original_value[2]
            new_matrix[1, 0] = original_value[3]
            new_matrix[1, 1] = original_value[4]
            new_matrix[1, 2] = original_value[5]
            self.assertMatrixAlmostEqual(self._GetValue(entity, Kratos.PK2_STRESS_TENSOR), new_matrix * 2, 12)

    def test_CornerCaseReshape(self) -> None:
        # Reshape a single scalar expression to itself
        input_expression = self._GetContainerExpression()
        self._Read(input_expression, Kratos.PRESSURE)
        combed = Kratos.Expression.Utils.Reshape(input_expression, [])
        array = Kratos.Vector(input_expression.GetExpression().NumberOfEntities() * input_expression.GetItemComponentCount())
        Kratos.Expression.CArrayExpressionIO.Write(combed, array)
        for i_entity, entity in enumerate(combed.GetContainer()):
            self.assertAlmostEqual(array[i_entity], self._GetValue(entity, Kratos.PRESSURE), 12)

    def test_Comb(self):
        a = self._GetContainerExpression()
        self._Read(a, Kratos.PRESSURE)
        b = self._GetContainerExpression()
        self._Read(b, Kratos.VELOCITY)

        combed = Kratos.Expression.Utils.Comb([a, b])
        combed *= 2.0
        self._Evaluate(combed, Kratos.PENALTY)

        for entity in a.GetContainer():
            original_p = self._GetValue(entity, Kratos.PRESSURE)
            original_v = self._GetValue(entity, Kratos.VELOCITY)
            new_vector = Kratos.Vector(4)
            new_vector[0] = original_p
            new_vector[1] = original_v[0]
            new_vector[2] = original_v[1]
            new_vector[3] = original_v[2]
            self.assertVectorAlmostEqual(self._GetValue(entity, Kratos.PENALTY), new_vector * 2, 12)

        combed = Kratos.Expression.Utils.Comb([a, b, a])
        combed *= 2.0
        self._Evaluate(combed, Kratos.PENALTY)

        for entity in a.GetContainer():
            original_p = self._GetValue(entity, Kratos.PRESSURE)
            original_v = self._GetValue(entity, Kratos.VELOCITY)
            new_vector = Kratos.Vector(5)
            new_vector[0] = original_p
            new_vector[1] = original_v[0]
            new_vector[2] = original_v[1]
            new_vector[3] = original_v[2]
            new_vector[4] = original_p
            self.assertVectorAlmostEqual(self._GetValue(entity, Kratos.PENALTY), new_vector * 2, 12)

    def test_CornerCaseComb(self) -> None:
        # Comb from a single scalar expression
        input_expression = self._GetContainerExpression()
        self._Read(input_expression, Kratos.PRESSURE)
        combed = Kratos.Expression.Utils.Comb([input_expression])
        array = Kratos.Vector(input_expression.GetExpression().NumberOfEntities() * input_expression.GetItemComponentCount())
        Kratos.Expression.CArrayExpressionIO.Write(combed, array)
        for i_entity, entity in enumerate(combed.GetContainer()):
            self.assertAlmostEqual(array[i_entity], self._GetValue(entity, Kratos.PRESSURE), 12)

    def test_SliceCombReshape(self):
        a = self._GetContainerExpression()
        self._Read(a, Kratos.INITIAL_STRAIN)

        km_comb = Kratos.Expression.Utils.Comb
        km_slice = Kratos.Expression.Utils.Slice
        km_reshape = Kratos.Expression.Utils.Reshape
        self._Evaluate(km_reshape(km_comb([a, km_slice(a, 2, 2), km_slice(a, 3, 2)]) * 2, [5, 2]), Kratos.PK2_STRESS_TENSOR)

        for entity in a.GetContainer():
            original_value = self._GetValue(entity, Kratos.INITIAL_STRAIN)
            new_matrix = Kratos.Matrix(5, 2)
            new_matrix[0, 0] = original_value[0]
            new_matrix[0, 1] = original_value[1]
            new_matrix[1, 0] = original_value[2]
            new_matrix[1, 1] = original_value[3]
            new_matrix[2, 0] = original_value[4]
            new_matrix[2, 1] = original_value[5]
            new_matrix[3, 0] = original_value[2]
            new_matrix[3, 1] = original_value[3]
            new_matrix[4, 0] = original_value[3]
            new_matrix[4, 1] = original_value[4]
            self.assertMatrixAlmostEqual(self._GetValue(entity, Kratos.PK2_STRESS_TENSOR), new_matrix * 2, 12)

    def test_GetContainer(self):
        a = self._GetContainerExpression()
        self.assertEqual(self._GetContainer(), a.GetContainer())

    def test_EmptyContainer(self):
        model_part = self.model.CreateModelPart("empty_model_part")
        a = type(self._GetContainerExpression())(model_part)
        self._Read(a, Kratos.VELOCITY)

        a *= 3
        self._Evaluate(a, Kratos.ACCELERATION)
        self.assertEqual(a.Evaluate().shape, (0, ))

    def test_VectorRead(self):
        a = self._GetContainerExpression()

        total_size = len(a.GetContainer()) * 3
        vector = Kratos.Vector(total_size)
        for i, entity in enumerate(a.GetContainer()):
            velocity = self._GetValue(entity, Kratos.VELOCITY)
            vector[i * 3] = velocity[0]
            vector[i * 3 + 1] = velocity[1]
            vector[i * 3 + 2] = velocity[2]

        Kratos.Expression.CArrayExpressionIO.Read(a, vector, [3])
        a *= 4.3
        self._Evaluate(a, Kratos.ACCELERATION)

        for entity in a.GetContainer():
            self.assertVectorAlmostEqual(self._GetValue(entity, Kratos.VELOCITY) * 4.3, self._GetValue(entity, Kratos.ACCELERATION))

    def test_VectorMove(self):
        a = self._GetContainerExpression()

        total_size = len(a.GetContainer()) * 3
        vector = Kratos.Vector(total_size)
        for i, entity in enumerate(a.GetContainer()):
            velocity = self._GetValue(entity, Kratos.VELOCITY)
            vector[i * 3] = velocity[0]
            vector[i * 3 + 1] = velocity[1]
            vector[i * 3 + 2] = velocity[2]

        Kratos.Expression.CArrayExpressionIO.Move(a, vector, [3])
        a *= 4.3
        self._Evaluate(a, Kratos.ACCELERATION)

        for i, entity in enumerate(a.GetContainer()):
            velocity = Kratos.Array3([vector[i * 3], vector[i * 3 + 1], vector[i * 3 + 2]])
            self.assertVectorAlmostEqual(velocity * 4.3, self._GetValue(entity, Kratos.ACCELERATION))

        for i in range(total_size):
            vector[i] *= 2

        self._Evaluate(a, Kratos.ACCELERATION)
        for i, entity in enumerate(a.GetContainer()):
            velocity = Kratos.Array3([vector[i * 3], vector[i * 3 + 1], vector[i * 3 + 2]])
            self.assertVectorAlmostEqual(velocity * 4.3, self._GetValue(entity, Kratos.ACCELERATION))

    def test_VectorWrite(self):
        a = self._GetContainerExpression()
        self._Read(a, Kratos.VELOCITY)

        vector = Kratos.Vector()
        Kratos.Expression.CArrayExpressionIO.Write(a, vector)

        for i, entity in enumerate(a.GetContainer()):
            velocity = Kratos.Array3([vector[i * 3], vector[i * 3 + 1], vector[i * 3 + 2]])
            self.assertVectorAlmostEqual(velocity, self._GetValue(entity, Kratos.VELOCITY))

    def test_GetMaxDepth(self):
        a = self._GetContainerExpression()
        self._Read(a, Kratos.GREEN_LAGRANGE_STRAIN_TENSOR)
        b = a + 10
        c = b * 2 + a
        d = c ** 2
        e = Kratos.Expression.Utils.Comb([a, a, d])
        f = Kratos.Expression.Utils.Reshape(e, [6, 2])
        g = Kratos.Expression.Utils.Slice(f, 2, 4)
        h = Kratos.Expression.Utils.Reshape(g, [2, 2])
        i = h - a
        j = Kratos.Expression.Utils.Abs(i)
        k = Kratos.Expression.Utils.EntitySum(j)
        self.assertEqual(k.GetMaxDepth(), 12)

    def test_Collapse(self):
        a = self._GetContainerExpression()
        self._Read(a, Kratos.GREEN_LAGRANGE_STRAIN_TENSOR)
        b = a + 10
        c = b * 2 + a
        d = c ** 2
        e = Kratos.Expression.Utils.Collapse(d)
        self.assertEqual(d.GetMaxDepth(), 5)
        self.assertEqual(e.GetMaxDepth(), 1)
        self.assertEqual(Kratos.Expression.Utils.NormInf(e-d), 0.0)

    def test_Abs(self):
        a = self._GetContainerExpression()
        self._Read(a, Kratos.PRESSURE)
        b = a * -1
        c = Kratos.Expression.Utils.Abs(b)
        for v1, v2, v3 in zip(a.Evaluate(), b.Evaluate(), c.Evaluate()):
            self.assertEqual(v1, -v2)
            self.assertEqual(v3, abs(v2))

    def test_EntityMin(self):
        a = self._GetContainerExpression()
        self._Read(a, Kratos.GREEN_LAGRANGE_STRAIN_TENSOR)
        b = Kratos.Expression.Utils.EntityMin(a)

        self.assertEqual(b.Evaluate().shape, (len(self._GetContainer()), ))
        for v1, v2 in zip(a.Evaluate(), b.Evaluate()):
            self.assertEqual(numpy.min(v1), v2)

    def test_EntityMax(self):
        a = self._GetContainerExpression()
        self._Read(a, Kratos.GREEN_LAGRANGE_STRAIN_TENSOR)
        b = Kratos.Expression.Utils.EntityMax(a * -1)

        self.assertEqual(b.Evaluate().shape, (len(self._GetContainer()), ))
        for v1, v2 in zip(a.Evaluate(), b.Evaluate()):
            self.assertEqual(numpy.max(-v1), v2)

    def test_EntitySum(self):
        a = self._GetContainerExpression()
        self._Read(a, Kratos.GREEN_LAGRANGE_STRAIN_TENSOR)
        b = Kratos.Expression.Utils.EntitySum(a)

        self.assertEqual(b.Evaluate().shape, (len(self._GetContainer()), ))
        for v1, v2 in zip(a.Evaluate(), b.Evaluate()):
            self.assertEqual(numpy.sum(v1), v2)

    def test_NormInf(self):
        a = self._GetContainerExpression()
        self._Read(a, Kratos.VELOCITY)
        a *= -1
        c = a.Evaluate().reshape([len(self._GetContainer()) * 3])
        if c.shape == (0,):
            # numpy norm inf throws an error if the array shape is zero.
            # hence this is checked before.
            norm_inf = 0.0
        else:
            norm_inf = numpy.linalg.norm(c, ord=numpy.inf)
        self.assertAlmostEqual(Kratos.Expression.Utils.NormInf(a), self.data_comm.MaxAll(norm_inf), 9)

    def test_NormL2(self):
        a = self._GetContainerExpression()
        self._Read(a, Kratos.VELOCITY)
        a *= -1
        c = a.Evaluate().reshape([len(self._GetContainer()) * 3])
        self.assertAlmostEqual(Kratos.Expression.Utils.NormL2(a), self.data_comm.SumAll(numpy.linalg.norm(c, ord=2) ** 2) ** 0.5, 9)

    def test_NormP(self):
        a = self._GetContainerExpression()
        self._Read(a, Kratos.VELOCITY)
        a *= -1
        c = a.Evaluate().reshape([len(self._GetContainer()) * 3])
        self.assertAlmostEqual(Kratos.Expression.Utils.NormP(a, 3), self.data_comm.SumAll(numpy.linalg.norm(c, ord=3) ** 3) ** (1/3), 9)

    def test_InnerProduct(self):
        a = self._GetContainerExpression()
        self._Read(a, Kratos.VELOCITY)
        self.assertAlmostEqual(Kratos.Expression.Utils.InnerProduct(a, a), self.data_comm.SumAll(numpy.linalg.norm(a.Evaluate()) ** 2), 9)

    @abstractmethod
    def _GetContainerExpression(self) -> ExpressionUnionType:
        pass

    @abstractmethod
    def _GetContainerType(self) -> Kratos.Globals.DataLocation:
        pass

    @abstractmethod
    def _GetContainer(self) -> Union[Kratos.NodesArray, Kratos.ConditionsArray, Kratos.ElementsArray]:
        pass

    @abstractmethod
    def _GetValue(self, entity, variable):
        pass

    @abstractmethod
    def _Read(self, container_expression, variable):
        pass

    @abstractmethod
    def _Evaluate(self, container_expression, variable):
        pass

class TestHistoricalContainerExpression(kratos_unittest.TestCase, TestContainerExpression):
    @classmethod
    def setUpClass(cls):
        cls.CreateEntities()

    def _GetContainerExpression(self) -> Kratos.Expression.NodalExpression:
        return Kratos.Expression.NodalExpression(self.model_part)

    def _GetContainerType(self) -> Kratos.Globals.DataLocation:
        return Kratos.Globals.DataLocation.NodeHistorical

    def _GetContainer(self):
        return self.model_part.GetCommunicator().LocalMesh().Nodes

    def _GetValue(self, entity, variable):
        return entity.GetSolutionStepValue(variable)

    def _Read(self, container_expression, variable):
        Kratos.Expression.VariableExpressionIO.Read(container_expression, variable, True)

    def _Evaluate(self, container_expression, variable):
        Kratos.Expression.VariableExpressionIO.Write(container_expression, variable, True)

class TestNodalContainerExpression(kratos_unittest.TestCase, TestContainerExpression):
    @classmethod
    def setUpClass(cls):
        cls.CreateEntities()

    def _GetContainerExpression(self) -> Kratos.Expression.NodalExpression:
        return Kratos.Expression.NodalExpression(self.model_part)

    def _GetContainerType(self) -> Kratos.Globals.DataLocation:
        return Kratos.Globals.DataLocation.NodeNonHistorical

    def _GetContainer(self):
        return self.model_part.GetCommunicator().LocalMesh().Nodes

    def _GetValue(self, entity, variable):
        return entity.GetValue(variable)

    def _Read(self, container_expression, variable):
        Kratos.Expression.VariableExpressionIO.Read(container_expression, variable, False)

    def _Evaluate(self, container_expression, variable):
        Kratos.Expression.VariableExpressionIO.Write(container_expression, variable, False)

class TestConditionContainerExpression(kratos_unittest.TestCase, TestContainerExpression):
    @classmethod
    def setUpClass(cls):
        cls.CreateEntities()

    def _GetContainerExpression(self) -> Kratos.Expression.ConditionExpression:
        return Kratos.Expression.ConditionExpression(self.model_part)

    def _GetContainerType(self) -> Kratos.Globals.DataLocation:
        return Kratos.Globals.DataLocation.Condition

    def _GetContainer(self):
        return self.model_part.GetCommunicator().LocalMesh().Conditions

    def _GetValue(self, entity, variable):
        return entity.GetValue(variable)

    def _Read(self, container_expression, variable):
        Kratos.Expression.VariableExpressionIO.Read(container_expression, variable)

    def _Evaluate(self, container_expression, variable):
        Kratos.Expression.VariableExpressionIO.Write(container_expression, variable)

    def testDomainSizeExpressionIOCondition(self):
        condition_exp = Kratos.Expression.ConditionExpression(self.model_part)
        Kratos.Expression.DomainSizeExpressionIO.Read(condition_exp)
        numpy_condition_exp = condition_exp.Evaluate()
        for i, condition in enumerate(self.model_part.Conditions):
            self.assertEqual(numpy_condition_exp[i], condition.GetGeometry().DomainSize())

    def testDomainSizeExpressionIOElement(self):
        element_exp = Kratos.Expression.ElementExpression(self.model_part)
        Kratos.Expression.DomainSizeExpressionIO.Read(element_exp)
        numpy_element_exp = element_exp.Evaluate()
        for i, element in enumerate(self.model_part.Elements):
            self.assertEqual(numpy_element_exp[i], element.GetGeometry().DomainSize())

    def testDomainSizeExpressionIOCondition_Empty(self):
        model = Kratos.Model()
        model_part = model.CreateModelPart("test")
        cond_exp = Kratos.Expression.ConditionExpression(model_part)
        Kratos.Expression.DomainSizeExpressionIO.Read(cond_exp)
        self.assertEqual(cond_exp.Evaluate().shape, (0,))

    def testDomainSizeExpressionIOElement_Empty(self):
        model = Kratos.Model()
        model_part = model.CreateModelPart("test")
        element_exp = Kratos.Expression.ElementExpression(model_part)
        Kratos.Expression.DomainSizeExpressionIO.Read(element_exp)
        self.assertEqual(element_exp.Evaluate().shape, (0,))

class TestElementContainerExpression(kratos_unittest.TestCase, TestContainerExpression):
    @classmethod
    def setUpClass(cls):
        cls.CreateEntities()

    def _GetContainerExpression(self):
        return Kratos.Expression.ElementExpression(self.model_part)

    def _GetContainerType(self) -> Kratos.Globals.DataLocation:
        return Kratos.Globals.DataLocation.Element

    def _GetContainer(self):
        return self.model_part.GetCommunicator().LocalMesh().Elements

    def _GetValue(self, entity, variable):
        return entity.GetValue(variable)

    def _Read(self, container_expression, variable):
        Kratos.Expression.VariableExpressionIO.Read(container_expression, variable)

    def _Evaluate(self, container_expression, variable):
        Kratos.Expression.VariableExpressionIO.Write(container_expression, variable)

class TestNodalPositionExpressionIO(kratos_unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.model = Kratos.Model()
        cls.model_part = cls.model.CreateModelPart("test")

        with kratos_unittest.WorkFolderScope(".", __file__, True):
            ReadModelPart("auxiliar_files_for_python_unittest/mdpa_files/two_dim_symmetrical_square", cls.model_part)

        for node in cls.model_part.Nodes:
            node.SetValue(Kratos.VELOCITY, Kratos.Array3([node.Id, node.Id + 1, node.Id + 2]))

    def test_WriteInitial(self):
        a = Kratos.Expression.NodalExpression(self.model_part)
        Kratos.Expression.VariableExpressionIO.Read(a, Kratos.VELOCITY, False)

        a *= 3

        Kratos.Expression.NodalPositionExpressionIO.Write(a, Kratos.Configuration.Initial)

        node: Kratos.Node
        for node in self.model_part.Nodes:
            self.assertVectorAlmostEqual(Kratos.Array3([node.X0, node.Y0, node.Z0]), node.GetValue(Kratos.VELOCITY) * 3)

    def test_WriteCurrent(self):
        a = Kratos.Expression.NodalExpression(self.model_part)
        Kratos.Expression.VariableExpressionIO.Read(a, Kratos.VELOCITY, False)

        a *= 4

        Kratos.Expression.NodalPositionExpressionIO.Write(a, Kratos.Configuration.Current)

        node: Kratos.Node
        for node in self.model_part.Nodes:
            self.assertVectorAlmostEqual(Kratos.Array3([node.X, node.Y, node.Z]), node.GetValue(Kratos.VELOCITY) * 4)

    def test_ReadInitial(self):
        a = Kratos.Expression.NodalExpression(self.model_part)
        Kratos.Expression.NodalPositionExpressionIO.Read(a, Kratos.Configuration.Initial)

        a *= 5
        Kratos.Expression.VariableExpressionIO.Write(a, Kratos.VELOCITY, False)


        node: Kratos.Node
        for node in self.model_part.Nodes:
            self.assertVectorAlmostEqual(Kratos.Array3([node.X0, node.Y0, node.Z0]) * 5, node.GetValue(Kratos.VELOCITY))

    def test_ReadCurrent(self):
        a = Kratos.Expression.NodalExpression(self.model_part)
        Kratos.Expression.NodalPositionExpressionIO.Read(a, Kratos.Configuration.Current)

        a *= 6
        Kratos.Expression.VariableExpressionIO.Write(a, Kratos.VELOCITY, False)

        node: Kratos.Node
        for node in self.model_part.Nodes:
            self.assertVectorAlmostEqual(Kratos.Array3([node.X0, node.Y0, node.Z0]) * 6, node.GetValue(Kratos.VELOCITY))


if __name__ == "__main__":
    Kratos.Tester.SetVerbosity(Kratos.Tester.Verbosity.PROGRESS)  # TESTS_OUTPUTS
    kratos_unittest.main()
