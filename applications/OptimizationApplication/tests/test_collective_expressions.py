import numpy
import KratosMultiphysics as Kratos
import KratosMultiphysics.OptimizationApplication as KratosOA
import KratosMultiphysics.KratosUnittest as kratos_unittest

class TestCollectiveExpressions(kratos_unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.model = Kratos.Model()
        cls.model_part = cls.model.CreateModelPart("test")
        cls.model_part.AddNodalSolutionStepVariable(Kratos.DENSITY)
        cls.model_part.AddNodalSolutionStepVariable(Kratos.PRESSURE)
        cls.model_part.AddNodalSolutionStepVariable(Kratos.ACCELERATION)
        cls.model_part.AddNodalSolutionStepVariable(Kratos.VELOCITY)

        number_of_nodes = 3
        for id in range(1, number_of_nodes + 1):
            node = cls.model_part.CreateNewNode(id, id, id+1, id+2)
            node.SetSolutionStepValue(Kratos.VELOCITY, Kratos.Array3([id+3, id+4, id+5]))
            node.SetSolutionStepValue(Kratos.PRESSURE, id+3)
            node.SetValue(Kratos.PRESSURE, id+3)
            node.SetValue(Kratos.VELOCITY, Kratos.Array3([id+3, id+4, id+5]))

        number_of_conditions = 4
        for id in range(1, number_of_conditions + 1):
            properties = cls.model_part.CreateNewProperties(id)
            properties.SetValue(Kratos.PRESSURE, id+400)
            properties.SetValue(Kratos.VELOCITY, Kratos.Array3([id+500, id+600, id+700]))
            condition = cls.model_part.CreateNewCondition("LineCondition2D2N", id + 1, [(id % number_of_nodes) + 1, ((id + 1) % number_of_nodes) + 1 ], properties)
            condition.SetValue(Kratos.PRESSURE, id+4)
            condition.SetValue(Kratos.VELOCITY, Kratos.Array3([id+5, id+6, id+7]))

        number_of_elements = 5
        for id in range(1, number_of_elements + 1):
            properties = cls.model_part.CreateNewProperties(id + number_of_conditions)
            properties.SetValue(Kratos.PRESSURE, id+500)
            properties.SetValue(Kratos.VELOCITY, Kratos.Array3([id+600, id+700, id+800]))
            element = cls.model_part.CreateNewElement("Element2D3N", id + 2, [(id % number_of_nodes) + 1, ((id + 1) % number_of_nodes) + 1, ((id + 2) % number_of_nodes) + 1 ], properties)
            element.SetValue(Kratos.PRESSURE, id+5)
            element.SetValue(Kratos.VELOCITY, Kratos.Array3([id+6, id+7, id+8]))

    def test_CollectiveExpressionsAdd(self):
        a = Kratos.Expression.HistoricalExpression(self.model_part)
        b = KratosOA.ContainerExpression.ElementPropertiesExpression(self.model_part)

        a.Read(Kratos.VELOCITY)
        b.Read(Kratos.PRESSURE)

        collective_1 = KratosOA.ContainerExpression.CollectiveExpressions()
        collective_1.Add(a)
        collective_1.Add(b)

        collective_2 = KratosOA.ContainerExpression.CollectiveExpressions()
        collective_2.Add(a)
        collective_2.Add(b)

        collective_3 = collective_1 * 10 + collective_2
        collective_3 += collective_1
        collective_4 = collective_3 + 4
        collective_4 += 10
        collective_4 += KratosOA.ContainerExpression.CollectiveExpressions(collective_1)

        collective_4.GetContainerExpressions()[0].Evaluate(Kratos.ACCELERATION)
        collective_4.GetContainerExpressions()[1].Evaluate(Kratos.DENSITY)

        for node in self.model_part.Nodes:
            self.assertVectorAlmostEqual(node.GetSolutionStepValue(Kratos.ACCELERATION), node.GetSolutionStepValue(Kratos.VELOCITY) * 13 + Kratos.Array3([14, 14, 14]), 12)
        for element in self.model_part.Elements:
            self.assertEqual(element.Properties[Kratos.DENSITY], element.Properties[Kratos.PRESSURE] * 13 + 14, 12)

    def test_CollectiveExpressionsSub(self):
        a = Kratos.Expression.HistoricalExpression(self.model_part)
        b = KratosOA.ContainerExpression.ElementPropertiesExpression(self.model_part)

        a.Read(Kratos.VELOCITY)
        b.Read(Kratos.PRESSURE)

        collective_1 = KratosOA.ContainerExpression.CollectiveExpressions()
        collective_1.Add(a)
        collective_1.Add(b)

        collective_2 = KratosOA.ContainerExpression.CollectiveExpressions()
        collective_2.Add(a)
        collective_2.Add(b)

        collective_3 = collective_1 * 10 - collective_2
        collective_3 -= collective_1
        collective_4 = collective_3 - 4
        collective_4 -= 10
        collective_4 -= KratosOA.ContainerExpression.CollectiveExpressions(collective_1)

        collective_4.GetContainerExpressions()[0].Evaluate(Kratos.ACCELERATION)
        collective_4.GetContainerExpressions()[1].Evaluate(Kratos.DENSITY)

        for node in self.model_part.Nodes:
            self.assertVectorAlmostEqual(node.GetSolutionStepValue(Kratos.ACCELERATION), node.GetSolutionStepValue(Kratos.VELOCITY) * 7 - Kratos.Array3([14, 14, 14]), 12)
        for element in self.model_part.Elements:
            self.assertEqual(element.Properties[Kratos.DENSITY], element.Properties[Kratos.PRESSURE] * 7 - 14, 12)

    def test_CollectiveExpressionsMul(self):
        a = Kratos.Expression.HistoricalExpression(self.model_part)
        b = KratosOA.ContainerExpression.ElementPropertiesExpression(self.model_part)

        a.Read(Kratos.VELOCITY)
        b.Read(Kratos.PRESSURE)

        collective_1 = KratosOA.ContainerExpression.CollectiveExpressions()
        collective_1.Add(a)
        collective_1.Add(b)

        collective_2 = KratosOA.ContainerExpression.CollectiveExpressions()
        collective_2.Add(a)
        collective_2.Add(b)

        collective_3 = collective_1 * 10
        collective_3 *= 2

        collective_3.GetContainerExpressions()[0].Evaluate(Kratos.ACCELERATION)
        collective_3.GetContainerExpressions()[1].Evaluate(Kratos.DENSITY)

        for node in self.model_part.Nodes:
            self.assertVectorAlmostEqual(node.GetSolutionStepValue(Kratos.ACCELERATION), node.GetSolutionStepValue(Kratos.VELOCITY) * 20, 12)
        for element in self.model_part.Elements:
            self.assertEqual(element.Properties[Kratos.DENSITY], element.Properties[Kratos.PRESSURE] * 20, 12)

    def test_CollectiveExpressionsDiv(self):
        a = Kratos.Expression.HistoricalExpression(self.model_part)
        b = KratosOA.ContainerExpression.ElementPropertiesExpression(self.model_part)

        a.Read(Kratos.VELOCITY)
        b.Read(Kratos.PRESSURE)

        collective_1 = KratosOA.ContainerExpression.CollectiveExpressions()
        collective_1.Add(a)
        collective_1.Add(b)

        collective_2 = KratosOA.ContainerExpression.CollectiveExpressions()
        collective_2.Add(a)
        collective_2.Add(b)

        collective_3 = collective_1 / 4
        collective_3 /= 2

        collective_3.GetContainerExpressions()[0].Evaluate(Kratos.ACCELERATION)
        collective_3.GetContainerExpressions()[1].Evaluate(Kratos.DENSITY)

        for node in self.model_part.Nodes:
            self.assertVectorAlmostEqual(node.GetSolutionStepValue(Kratos.ACCELERATION), node.GetSolutionStepValue(Kratos.VELOCITY) / 8, 12)
        for element in self.model_part.Elements:
            self.assertEqual(element.Properties[Kratos.DENSITY], element.Properties[Kratos.PRESSURE] / 8, 12)

    def test_CollectiveExpressionsPow(self):
        a = Kratos.Expression.HistoricalExpression(self.model_part)
        b = KratosOA.ContainerExpression.ElementPropertiesExpression(self.model_part)

        a.Read(Kratos.VELOCITY)
        b.Read(Kratos.PRESSURE)

        collective_1 = KratosOA.ContainerExpression.CollectiveExpressions()
        collective_1.Add(a)
        collective_1.Add(b)

        collective_2 = KratosOA.ContainerExpression.CollectiveExpressions()
        collective_2.Add(a)
        collective_2.Add(b)

        collective_3 = collective_1 ** (collective_1 / 1e+3)
        collective_3 **= 2

        collective_3.GetContainerExpressions()[0].Evaluate(Kratos.ACCELERATION)
        collective_3.GetContainerExpressions()[1].Evaluate(Kratos.DENSITY)

        for node in self.model_part.Nodes:
            v = node.GetSolutionStepValue(Kratos.VELOCITY)
            self.assertVectorAlmostEqual(node.GetSolutionStepValue(Kratos.ACCELERATION), Kratos.Array3([v[0]**(2*v[0]/1e+3), v[1]**(2*v[1]/1e+3), v[2]**(2*v[2]/1e+3)]) , 12)
        for element in self.model_part.Elements:
            self.assertAlmostEqual(element.Properties[Kratos.DENSITY], element.Properties[Kratos.PRESSURE] ** (2 * element.Properties[Kratos.PRESSURE] / 1e+3), 12)

    def test_CollectiveExpressionsNeg(self):
        a = Kratos.Expression.HistoricalExpression(self.model_part)
        b = KratosOA.ContainerExpression.ElementPropertiesExpression(self.model_part)

        a.Read(Kratos.VELOCITY)
        b.Read(Kratos.PRESSURE)

        collective_1 = KratosOA.ContainerExpression.CollectiveExpressions()
        collective_1.Add(a)
        collective_1.Add(b)

        collective_2 = KratosOA.ContainerExpression.CollectiveExpressions()
        collective_2.Add(a)
        collective_2.Add(b)

        collective_3 = -collective_1

        collective_3.GetContainerExpressions()[0].Evaluate(Kratos.ACCELERATION)
        collective_3.GetContainerExpressions()[1].Evaluate(Kratos.DENSITY)

        for node in self.model_part.Nodes:
            v = node.GetSolutionStepValue(Kratos.VELOCITY)
            self.assertVectorAlmostEqual(node.GetSolutionStepValue(Kratos.ACCELERATION), v*(-1) , 12)
        for element in self.model_part.Elements:
            self.assertEqual(element.Properties[Kratos.DENSITY], -element.Properties[Kratos.PRESSURE], 12)

    def test_IsCompatibleWith(self):
        additional_model_part = self.model.CreateModelPart("additional_model_part")

        diff_size_model_part = self.model.CreateModelPart("diff_size_model_part")

        for properties in self.model_part.GetProperties():
            additional_model_part.AddProperties(properties)

        for node in self.model_part.Nodes:
            additional_model_part.AddNode(node)

        for condition in self.model_part.Conditions:
            additional_model_part.AddCondition(condition)

        for element in self.model_part.Elements:
            additional_model_part.AddElement(element)

        a = Kratos.Expression.HistoricalExpression(self.model_part)
        b = KratosOA.ContainerExpression.ElementPropertiesExpression(self.model_part)
        c = Kratos.Expression.HistoricalExpression(additional_model_part)
        d = Kratos.Expression.HistoricalExpression(diff_size_model_part)

        a.Read(Kratos.VELOCITY)
        b.Read(Kratos.PRESSURE)
        c.Read(Kratos.VELOCITY)
        d.Read(Kratos.VELOCITY)

        collective_1 = KratosOA.ContainerExpression.CollectiveExpressions([a, b])
        self.assertTrue(collective_1.IsCompatibleWith(KratosOA.ContainerExpression.CollectiveExpressions([a, b])))
        self.assertFalse(collective_1.IsCompatibleWith(KratosOA.ContainerExpression.CollectiveExpressions([a])))
        self.assertFalse(collective_1.IsCompatibleWith(KratosOA.ContainerExpression.CollectiveExpressions([b, a])))
        self.assertTrue(collective_1.IsCompatibleWith(KratosOA.ContainerExpression.CollectiveExpressions([c, b])))
        self.assertFalse(collective_1.IsCompatibleWith(KratosOA.ContainerExpression.CollectiveExpressions([c, d])))

    def test_ReadEvaluate1(self):
        a = Kratos.Expression.HistoricalExpression(self.model_part)
        b = Kratos.Expression.NodalNonHistoricalExpression(self.model_part)
        c = Kratos.Expression.ConditionNonHistoricalExpression(self.model_part)
        d = Kratos.Expression.ElementNonHistoricalExpression(self.model_part)

        collective = KratosOA.ContainerExpression.CollectiveExpressions([a, b, c, d])

        # reads VELOCITY from different container expressions
        collective.Read(Kratos.VELOCITY)

        # do some calculations on all container expressions. these are lazy expressions, hence light wieght operations.
        collective += 2
        collective *= 3

        # finally evaluate the lazy expressions and put it to a numpy continuous array
        result = collective.Evaluate()
        self.assertEqual(result.shape, (self.model_part.NumberOfNodes() * 2 * 3 + self.model_part.NumberOfConditions() * 3 + self.model_part.NumberOfElements() * 3, ))

        # now check
        index = 0
        for node in self.model_part.Nodes:
            velocity = node.GetSolutionStepValue(Kratos.VELOCITY)
            self.assertEqual((velocity[0] + 2) * 3, result[index])
            self.assertEqual((velocity[1] + 2) * 3, result[index + 1])
            self.assertEqual((velocity[2] + 2) * 3, result[index + 2])
            index += 3

        for node in self.model_part.Nodes:
            velocity = node.GetValue(Kratos.VELOCITY)
            self.assertEqual((velocity[0] + 2) * 3, result[index])
            self.assertEqual((velocity[1] + 2) * 3, result[index + 1])
            self.assertEqual((velocity[2] + 2) * 3, result[index + 2])
            index += 3

        for condition in self.model_part.Conditions:
            velocity = condition.GetValue(Kratos.VELOCITY)
            self.assertEqual((velocity[0] + 2) * 3, result[index])
            self.assertEqual((velocity[1] + 2) * 3, result[index + 1])
            self.assertEqual((velocity[2] + 2) * 3, result[index + 2])
            index += 3

        for element in self.model_part.Elements:
            velocity = element.GetValue(Kratos.VELOCITY)
            self.assertEqual((velocity[0] + 2) * 3, result[index])
            self.assertEqual((velocity[1] + 2) * 3, result[index + 1])
            self.assertEqual((velocity[2] + 2) * 3, result[index + 2])
            index += 3

    def test_ReadEvaluate2(self):
        a = Kratos.Expression.HistoricalExpression(self.model_part)
        b = Kratos.Expression.NodalNonHistoricalExpression(self.model_part)
        c = Kratos.Expression.ConditionNonHistoricalExpression(self.model_part)
        d = Kratos.Expression.ElementNonHistoricalExpression(self.model_part)

        collective = KratosOA.ContainerExpression.CollectiveExpressions([a, b, c, d])

        # reads VELOCITY from different container expressions
        collective.Read([Kratos.VELOCITY, Kratos.PRESSURE, Kratos.VELOCITY, Kratos.PRESSURE])

        # do some calculations on all container expressions. these are lazy expressions, hence light wieght operations.
        collective += 2
        collective *= 3

        # finally evaluate the lazy expressions and put it to a numpy continuous array
        result = collective.Evaluate()
        self.assertEqual(result.shape, (self.model_part.NumberOfNodes() * 4 + self.model_part.NumberOfConditions() * 3 + self.model_part.NumberOfElements(), ))

        # now check
        index = 0
        for node in self.model_part.Nodes:
            velocity = node.GetSolutionStepValue(Kratos.VELOCITY)
            self.assertEqual((velocity[0] + 2) * 3, result[index])
            self.assertEqual((velocity[1] + 2) * 3, result[index + 1])
            self.assertEqual((velocity[2] + 2) * 3, result[index + 2])
            index += 3

        for node in self.model_part.Nodes:
            pressure = node.GetValue(Kratos.PRESSURE)
            self.assertEqual((pressure + 2) * 3, result[index])
            index += 1

        for condition in self.model_part.Conditions:
            velocity = condition.GetValue(Kratos.VELOCITY)
            self.assertEqual((velocity[0] + 2) * 3, result[index])
            self.assertEqual((velocity[1] + 2) * 3, result[index + 1])
            self.assertEqual((velocity[2] + 2) * 3, result[index + 2])
            index += 3

        for element in self.model_part.Elements:
            pressure = element.GetValue(Kratos.PRESSURE)
            self.assertEqual((pressure + 2) * 3, result[index])
            index += 1

    def test_ReadEvaluate3(self):
        a = Kratos.Expression.HistoricalExpression(self.model_part)
        b = Kratos.Expression.NodalNonHistoricalExpression(self.model_part)
        c = Kratos.Expression.ConditionNonHistoricalExpression(self.model_part)
        d = Kratos.Expression.ElementNonHistoricalExpression(self.model_part)

        collective = KratosOA.ContainerExpression.CollectiveExpressions([a, b, c, d])

        number_of_entities = numpy.array([self.model_part.NumberOfNodes(), self.model_part.NumberOfNodes(), self.model_part.NumberOfConditions(), self.model_part.NumberOfElements()])
        data = numpy.arange(1, self.model_part.NumberOfNodes() * 4 + self.model_part.NumberOfConditions() * 3 + self.model_part.NumberOfElements() + 1, 1.0)

        # here we copy data
        collective.Read(data, number_of_entities, [[3], [], [3], []])

        # do some calculations on all container expressions. these are lazy expressions, hence light wieght operations.
        collective += 2
        collective *= 3

        collective.Evaluate([Kratos.ACCELERATION, Kratos.DENSITY, Kratos.ACCELERATION, Kratos.DENSITY])

        # now check
        index = 0
        for node in self.model_part.Nodes:
            acceleration = node.GetSolutionStepValue(Kratos.ACCELERATION)
            self.assertEqual(acceleration[0], (data[index] + 2) * 3)
            self.assertEqual(acceleration[1], (data[index + 1] + 2) * 3)
            self.assertEqual(acceleration[2], (data[index + 2] + 2) * 3)
            index += 3

        for node in self.model_part.Nodes:
            pressure = node.GetValue(Kratos.DENSITY)
            self.assertEqual(pressure, (data[index] + 2) * 3)
            index += 1

        for condition in self.model_part.Conditions:
            acceleration = condition.GetValue(Kratos.ACCELERATION)
            self.assertEqual(acceleration[0], (data[index] + 2) * 3)
            self.assertEqual(acceleration[1], (data[index + 1] + 2) * 3)
            self.assertEqual(acceleration[2], (data[index + 2] + 2) * 3)
            index += 3

        for element in self.model_part.Elements:
            pressure = element.GetValue(Kratos.DENSITY)
            self.assertEqual(pressure, (data[index] + 2) * 3)
            index += 1

    def test_Move(self):
        a = Kratos.Expression.HistoricalExpression(self.model_part)
        b = Kratos.Expression.NodalNonHistoricalExpression(self.model_part)
        c = Kratos.Expression.ConditionNonHistoricalExpression(self.model_part)
        d = Kratos.Expression.ElementNonHistoricalExpression(self.model_part)

        collective = KratosOA.ContainerExpression.CollectiveExpressions([a, b, c, d])

        number_of_entities = [self.model_part.NumberOfNodes(), self.model_part.NumberOfNodes(), self.model_part.NumberOfConditions(), self.model_part.NumberOfElements()]
        data = numpy.arange(1, self.model_part.NumberOfNodes() * 4 + self.model_part.NumberOfConditions() * 3 + self.model_part.NumberOfElements() + 1, 1.0)

        # here we move data. The life time is not managed by the collective. If moved data is destroyed, then use of collective can seg fault.
        collective.MoveFrom(data, number_of_entities, [[3], [], [3], []])

        # do some calculations on all container expressions. these are lazy expressions, hence light wieght operations.
        collective += 2
        collective *= 3

        collective.Evaluate([Kratos.ACCELERATION, Kratos.DENSITY, Kratos.ACCELERATION, Kratos.DENSITY])

        # now check for initial values
        index = 0
        for node in self.model_part.Nodes:
            acceleration = node.GetSolutionStepValue(Kratos.ACCELERATION)
            self.assertEqual(acceleration[0], (data[index] + 2) * 3)
            self.assertEqual(acceleration[1], (data[index + 1] + 2) * 3)
            self.assertEqual(acceleration[2], (data[index + 2] + 2) * 3)
            index += 3

        for node in self.model_part.Nodes:
            pressure = node.GetValue(Kratos.DENSITY)
            self.assertEqual(pressure, (data[index] + 2) * 3)
            index += 1

        for condition in self.model_part.Conditions:
            acceleration = condition.GetValue(Kratos.ACCELERATION)
            self.assertEqual(acceleration[0], (data[index] + 2) * 3)
            self.assertEqual(acceleration[1], (data[index + 1] + 2) * 3)
            self.assertEqual(acceleration[2], (data[index + 2] + 2) * 3)
            index += 3

        for element in self.model_part.Elements:
            pressure = element.GetValue(Kratos.DENSITY)
            self.assertEqual(pressure, (data[index] + 2) * 3)
            index += 1

        # now change the numpy array to check whether we are still referening to data from numpy arraydata
        data += 1.0
        collective.Evaluate([Kratos.ACCELERATION, Kratos.DENSITY, Kratos.ACCELERATION, Kratos.DENSITY])

        # now check for changed values
        index = 0
        for node in self.model_part.Nodes:
            acceleration = node.GetSolutionStepValue(Kratos.ACCELERATION)
            self.assertEqual(acceleration[0], (data[index] + 2) * 3)
            self.assertEqual(acceleration[1], (data[index + 1] + 2) * 3)
            self.assertEqual(acceleration[2], (data[index + 2] + 2) * 3)
            index += 3

        for node in self.model_part.Nodes:
            pressure = node.GetValue(Kratos.DENSITY)
            self.assertEqual(pressure, (data[index] + 2) * 3)
            index += 1

        for condition in self.model_part.Conditions:
            acceleration = condition.GetValue(Kratos.ACCELERATION)
            self.assertEqual(acceleration[0], (data[index] + 2) * 3)
            self.assertEqual(acceleration[1], (data[index + 1] + 2) * 3)
            self.assertEqual(acceleration[2], (data[index + 2] + 2) * 3)
            index += 3

        for element in self.model_part.Elements:
            pressure = element.GetValue(Kratos.DENSITY)
            self.assertEqual(pressure, (data[index] + 2) * 3)
            index += 1

    def test_ReadMoveErrors(self):
        a = Kratos.Expression.HistoricalExpression(self.model_part)
        b = Kratos.Expression.NodalNonHistoricalExpression(self.model_part)
        c = Kratos.Expression.ConditionNonHistoricalExpression(self.model_part)
        d = Kratos.Expression.ElementNonHistoricalExpression(self.model_part)

        collective = KratosOA.ContainerExpression.CollectiveExpressions([a, b, c, d])

        number_of_entities = numpy.array([self.model_part.NumberOfNodes(), self.model_part.NumberOfNodes(), self.model_part.NumberOfConditions(), self.model_part.NumberOfElements()])
        total_entities = self.model_part.NumberOfNodes() * 4 + self.model_part.NumberOfConditions() * 3 + self.model_part.NumberOfElements() + 1

        with self.assertRaises(TypeError):
            numpy_array = numpy.arange(0, total_entities)
            collective.Read(numpy_array, number_of_entities, [[3], [], [3], []])

        with self.assertRaises(TypeError):
            numpy_array = numpy.arange(0, total_entities, dtype=numpy.float32)
            collective.Read(numpy_array, number_of_entities, [[3], [], [3], []])

        with self.assertRaises(TypeError):
            numpy_array = numpy.arange(0, total_entities, dtype=numpy.int32)
            collective.Read(numpy_array, number_of_entities, [[3], [], [3], []])

        with self.assertRaises(TypeError):
            numpy_array = numpy.arange(0, total_entities, dtype=numpy.int64)
            collective.Read(numpy_array, number_of_entities, [[3], [], [3], []])

        with self.assertRaises(TypeError):
            numpy_array = numpy.arange(0, total_entities)
            collective.MoveFrom(numpy_array, number_of_entities, [[3], [], [3], []])

        with self.assertRaises(TypeError):
            numpy_array = numpy.arange(0, total_entities, dtype=numpy.float32)
            collective.MoveFrom(numpy_array, number_of_entities, [[3], [], [3], []])

        with self.assertRaises(TypeError):
            numpy_array = numpy.arange(0, total_entities, dtype=numpy.int32)
            collective.MoveFrom(numpy_array, number_of_entities, [[3], [], [3], []])

        with self.assertRaises(TypeError):
            numpy_array = numpy.arange(0, total_entities, dtype=numpy.int64)
            collective.MoveFrom(numpy_array, number_of_entities, [[3], [], [3], []])

if __name__ == "__main__":
    Kratos.Tester.SetVerbosity(Kratos.Tester.Verbosity.PROGRESS)  # TESTS_OUTPUTS
    kratos_unittest.main()