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
        a = Kratos.Expression.NodalExpression(self.model_part)
        b = Kratos.Expression.ElementExpression(self.model_part)

        Kratos.Expression.VariableExpressionIO.Read(a, Kratos.VELOCITY, True)
        KratosOA.PropertiesVariableExpressionIO.Read(b, Kratos.PRESSURE)

        collective_1 = KratosOA.CollectiveExpression()
        collective_1.Add(a)
        collective_1.Add(b)

        collective_2 = KratosOA.CollectiveExpression()
        collective_2.Add(a)
        collective_2.Add(b)

        collective_3 = collective_1 * 10 + collective_2
        collective_3 += collective_1
        collective_4 = collective_3 + 4
        collective_4 += 10
        collective_4 += collective_1.Clone()

        Kratos.Expression.VariableExpressionIO.Write(collective_4.GetContainerExpressions()[0], Kratos.ACCELERATION, True)
        KratosOA.PropertiesVariableExpressionIO.Write(collective_4.GetContainerExpressions()[1], Kratos.DENSITY)

        for node in self.model_part.Nodes:
            self.assertVectorAlmostEqual(node.GetSolutionStepValue(Kratos.ACCELERATION), node.GetSolutionStepValue(Kratos.VELOCITY) * 13 + Kratos.Array3([14, 14, 14]), 12)
        for element in self.model_part.Elements:
            self.assertEqual(element.Properties[Kratos.DENSITY], element.Properties[Kratos.PRESSURE] * 13 + 14, 12)

    def test_CollectiveExpressionsSub(self):
        a = Kratos.Expression.NodalExpression(self.model_part)
        b = Kratos.Expression.ElementExpression(self.model_part)

        Kratos.Expression.VariableExpressionIO.Read(a, Kratos.VELOCITY, True)
        KratosOA.PropertiesVariableExpressionIO.Read(b, Kratos.PRESSURE)

        collective_1 = KratosOA.CollectiveExpression()
        collective_1.Add(a)
        collective_1.Add(b)

        collective_2 = KratosOA.CollectiveExpression()
        collective_2.Add(a)
        collective_2.Add(b)

        collective_3 = collective_1 * 10 - collective_2
        collective_3 -= collective_1
        collective_4 = collective_3 - 4
        collective_4 -= 10
        collective_4 -= collective_1.Clone()

        Kratos.Expression.VariableExpressionIO.Write(collective_4.GetContainerExpressions()[0], Kratos.ACCELERATION, True)
        KratosOA.PropertiesVariableExpressionIO.Write(collective_4.GetContainerExpressions()[1], Kratos.DENSITY)

        for node in self.model_part.Nodes:
            self.assertVectorAlmostEqual(node.GetSolutionStepValue(Kratos.ACCELERATION), node.GetSolutionStepValue(Kratos.VELOCITY) * 7 - Kratos.Array3([14, 14, 14]), 12)
        for element in self.model_part.Elements:
            self.assertEqual(element.Properties[Kratos.DENSITY], element.Properties[Kratos.PRESSURE] * 7 - 14, 12)

    def test_CollectiveExpressionsMul(self):
        a = Kratos.Expression.NodalExpression(self.model_part)
        b = Kratos.Expression.ElementExpression(self.model_part)

        Kratos.Expression.VariableExpressionIO.Read(a, Kratos.VELOCITY, True)
        KratosOA.PropertiesVariableExpressionIO.Read(b, Kratos.PRESSURE)

        collective_1 = KratosOA.CollectiveExpression()
        collective_1.Add(a)
        collective_1.Add(b)

        collective_2 = KratosOA.CollectiveExpression()
        collective_2.Add(a)
        collective_2.Add(b)

        collective_3 = collective_1 * 10
        collective_3 *= 2

        Kratos.Expression.VariableExpressionIO.Write(collective_3.GetContainerExpressions()[0], Kratos.ACCELERATION, True)
        KratosOA.PropertiesVariableExpressionIO.Write(collective_3.GetContainerExpressions()[1], Kratos.DENSITY)

        for node in self.model_part.Nodes:
            self.assertVectorAlmostEqual(node.GetSolutionStepValue(Kratos.ACCELERATION), node.GetSolutionStepValue(Kratos.VELOCITY) * 20, 12)
        for element in self.model_part.Elements:
            self.assertEqual(element.Properties[Kratos.DENSITY], element.Properties[Kratos.PRESSURE] * 20, 12)

    def test_CollectiveExpressionsDiv(self):
        a = Kratos.Expression.NodalExpression(self.model_part)
        b = Kratos.Expression.ElementExpression(self.model_part)

        Kratos.Expression.VariableExpressionIO.Read(a, Kratos.VELOCITY, True)
        KratosOA.PropertiesVariableExpressionIO.Read(b, Kratos.PRESSURE)

        collective_1 = KratosOA.CollectiveExpression()
        collective_1.Add(a)
        collective_1.Add(b)

        collective_2 = KratosOA.CollectiveExpression()
        collective_2.Add(a)
        collective_2.Add(b)

        collective_3 = collective_1 / 4
        collective_3 /= 2

        Kratos.Expression.VariableExpressionIO.Write(collective_3.GetContainerExpressions()[0], Kratos.ACCELERATION, True)
        KratosOA.PropertiesVariableExpressionIO.Write(collective_3.GetContainerExpressions()[1], Kratos.DENSITY)

        for node in self.model_part.Nodes:
            self.assertVectorAlmostEqual(node.GetSolutionStepValue(Kratos.ACCELERATION), node.GetSolutionStepValue(Kratos.VELOCITY) / 8, 12)
        for element in self.model_part.Elements:
            self.assertEqual(element.Properties[Kratos.DENSITY], element.Properties[Kratos.PRESSURE] / 8, 12)

    def test_CollectiveExpressionsPow(self):
        a = Kratos.Expression.NodalExpression(self.model_part)
        b = Kratos.Expression.ElementExpression(self.model_part)

        Kratos.Expression.VariableExpressionIO.Read(a, Kratos.VELOCITY, True)
        KratosOA.PropertiesVariableExpressionIO.Read(b, Kratos.PRESSURE)

        collective_1 = KratosOA.CollectiveExpression()
        collective_1.Add(a)
        collective_1.Add(b)

        collective_2 = KratosOA.CollectiveExpression()
        collective_2.Add(a)
        collective_2.Add(b)

        collective_3 = collective_1 ** (collective_1 / 1e+3)
        collective_3 **= 2

        Kratos.Expression.VariableExpressionIO.Write(collective_3.GetContainerExpressions()[0], Kratos.ACCELERATION, True)
        KratosOA.PropertiesVariableExpressionIO.Write(collective_3.GetContainerExpressions()[1], Kratos.DENSITY)

        for node in self.model_part.Nodes:
            v = node.GetSolutionStepValue(Kratos.VELOCITY)
            self.assertVectorAlmostEqual(node.GetSolutionStepValue(Kratos.ACCELERATION), Kratos.Array3([v[0]**(2*v[0]/1e+3), v[1]**(2*v[1]/1e+3), v[2]**(2*v[2]/1e+3)]) , 12)
        for element in self.model_part.Elements:
            self.assertAlmostEqual(element.Properties[Kratos.DENSITY], element.Properties[Kratos.PRESSURE] ** (2 * element.Properties[Kratos.PRESSURE] / 1e+3), 12)

    def test_CollectiveExpressionsNeg(self):
        a = Kratos.Expression.NodalExpression(self.model_part)
        b = Kratos.Expression.ElementExpression(self.model_part)

        Kratos.Expression.VariableExpressionIO.Read(a, Kratos.VELOCITY, True)
        KratosOA.PropertiesVariableExpressionIO.Read(b, Kratos.PRESSURE)

        collective_1 = KratosOA.CollectiveExpression()
        collective_1.Add(a)
        collective_1.Add(b)

        collective_2 = KratosOA.CollectiveExpression()
        collective_2.Add(a)
        collective_2.Add(b)

        collective_3 = -collective_1

        Kratos.Expression.VariableExpressionIO.Write(collective_3.GetContainerExpressions()[0], Kratos.ACCELERATION, True)
        KratosOA.PropertiesVariableExpressionIO.Write(collective_3.GetContainerExpressions()[1], Kratos.DENSITY)

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

        a = Kratos.Expression.NodalExpression(self.model_part)
        b = Kratos.Expression.ElementExpression(self.model_part)
        c = Kratos.Expression.NodalExpression(additional_model_part)
        d = Kratos.Expression.NodalExpression(diff_size_model_part)

        Kratos.Expression.VariableExpressionIO.Read(a, Kratos.VELOCITY, True)
        Kratos.Expression.VariableExpressionIO.Read(b, Kratos.PRESSURE)
        Kratos.Expression.VariableExpressionIO.Read(c, Kratos.VELOCITY, True)
        Kratos.Expression.VariableExpressionIO.Read(d, Kratos.VELOCITY, True)

        collective_1 = KratosOA.CollectiveExpression([a, b])
        self.assertTrue(collective_1.IsCompatibleWith(KratosOA.CollectiveExpression([a, b])))
        self.assertFalse(collective_1.IsCompatibleWith(KratosOA.CollectiveExpression([a])))
        self.assertFalse(collective_1.IsCompatibleWith(KratosOA.CollectiveExpression([b, a])))
        self.assertTrue(collective_1.IsCompatibleWith(KratosOA.CollectiveExpression([c, b])))
        self.assertFalse(collective_1.IsCompatibleWith(KratosOA.CollectiveExpression([c, d])))

    def test_ReadEvaluate1(self):
        a = Kratos.Expression.NodalExpression(self.model_part)
        b = Kratos.Expression.NodalExpression(self.model_part)
        c = Kratos.Expression.ConditionExpression(self.model_part)
        d = Kratos.Expression.ElementExpression(self.model_part)

        collective = KratosOA.CollectiveExpression([a, b, c, d])

        # reads VELOCITY from different container expressions
        ceio = KratosOA.CollectiveExpressionIO
        historical = ceio.HistoricalVariable(Kratos.VELOCITY)
        non_historical = ceio.NonHistoricalVariable(Kratos.VELOCITY)
        KratosOA.CollectiveExpressionIO.Read(collective,  [historical, non_historical, non_historical, non_historical])

        # do some calculations on all container expressions. these are lazy expressions, hence light wieght operations.
        collective += 2
        collective *= 3

        # finally evaluate the lazy expressions and put it to a numpy continuous array
        result = collective.Evaluate()
        self.assertEqual(result.shape, (self.model_part.NumberOfNodes() * 2 + self.model_part.NumberOfConditions() + self.model_part.NumberOfElements(), 3))

        # now check
        offset = 0
        for i, node in enumerate(self.model_part.Nodes):
            velocity = node.GetSolutionStepValue(Kratos.VELOCITY)
            self.assertEqual((velocity[0] + 2) * 3, result[i + offset, 0])
            self.assertEqual((velocity[1] + 2) * 3, result[i + offset, 1])
            self.assertEqual((velocity[2] + 2) * 3, result[i + offset, 2])

        offset += self.model_part.NumberOfNodes()
        for i, node in enumerate(self.model_part.Nodes):
            velocity = node.GetValue(Kratos.VELOCITY)
            self.assertEqual((velocity[0] + 2) * 3, result[i + offset, 0])
            self.assertEqual((velocity[1] + 2) * 3, result[i + offset, 1])
            self.assertEqual((velocity[2] + 2) * 3, result[i + offset, 2])

        offset += self.model_part.NumberOfNodes()
        for i, condition in enumerate(self.model_part.Conditions):
            velocity = condition.GetValue(Kratos.VELOCITY)
            self.assertEqual((velocity[0] + 2) * 3, result[i + offset, 0])
            self.assertEqual((velocity[1] + 2) * 3, result[i + offset, 1])
            self.assertEqual((velocity[2] + 2) * 3, result[i + offset, 2])

        offset += self.model_part.NumberOfConditions()
        for i, element in enumerate(self.model_part.Elements):
            velocity = element.GetValue(Kratos.VELOCITY)
            self.assertEqual((velocity[0] + 2) * 3, result[i + offset, 0])
            self.assertEqual((velocity[1] + 2) * 3, result[i + offset, 1])
            self.assertEqual((velocity[2] + 2) * 3, result[i + offset, 2])

    def test_ReadEvaluate2(self):
        a = Kratos.Expression.NodalExpression(self.model_part)
        b = Kratos.Expression.NodalExpression(self.model_part)
        c = Kratos.Expression.ConditionExpression(self.model_part)
        d = Kratos.Expression.ElementExpression(self.model_part)

        collective = KratosOA.CollectiveExpression([a, b, c, d])

        # reads VELOCITY from different container expressions
        ceio = KratosOA.CollectiveExpressionIO
        historical = ceio.HistoricalVariable
        non_historical = ceio.NonHistoricalVariable
        KratosOA.CollectiveExpressionIO.Read(collective, [historical(Kratos.VELOCITY), non_historical(Kratos.PRESSURE), non_historical(Kratos.VELOCITY), non_historical(Kratos.PRESSURE)])

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
        a = Kratos.Expression.NodalExpression(self.model_part)
        b = Kratos.Expression.NodalExpression(self.model_part)
        c = Kratos.Expression.ConditionExpression(self.model_part)
        d = Kratos.Expression.ElementExpression(self.model_part)

        collective = KratosOA.CollectiveExpression([a, b, c, d])

        data = numpy.arange(1, self.model_part.NumberOfNodes() * 4 + self.model_part.NumberOfConditions() * 3 + self.model_part.NumberOfElements() + 1, 1.0)

        # here we copy data
        KratosOA.CollectiveExpressionIO.Read(collective, data, [[3], [], [3], []])

        # do some calculations on all container expressions. these are lazy expressions, hence light wieght operations.
        collective += 2
        collective *= 3

        ceio = KratosOA.CollectiveExpressionIO
        historical = ceio.HistoricalVariable
        non_historical = ceio.NonHistoricalVariable
        KratosOA.CollectiveExpressionIO.Write(collective, [historical(Kratos.ACCELERATION), non_historical(Kratos.DENSITY), non_historical(Kratos.ACCELERATION), non_historical(Kratos.DENSITY)])

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

    def test_ReadEvaluate4(self):
        a = Kratos.Expression.NodalExpression(self.model_part)
        b = Kratos.Expression.NodalExpression(self.model_part)
        c = Kratos.Expression.ConditionExpression(self.model_part)
        d = Kratos.Expression.ElementExpression(self.model_part)

        collective = KratosOA.CollectiveExpression([a, b, c, d])

        data = numpy.arange(1, self.model_part.NumberOfNodes() * 6 + self.model_part.NumberOfConditions() * 3 + self.model_part.NumberOfElements() * 3 + 1, 1.0)
        data = data.reshape((self.model_part.NumberOfNodes() * 2 + self.model_part.NumberOfConditions() + self.model_part.NumberOfElements(), 3))

        # here we copy data
        KratosOA.CollectiveExpressionIO.Read(collective, data)

        # do some calculations on all container expressions. these are lazy expressions, hence light wieght operations.
        collective += 2
        collective *= 3

        ceio = KratosOA.CollectiveExpressionIO
        historical = ceio.HistoricalVariable
        non_historical = ceio.NonHistoricalVariable
        KratosOA.CollectiveExpressionIO.Write(collective, [historical(Kratos.ACCELERATION), non_historical(Kratos.ACCELERATION), non_historical(Kratos.ACCELERATION), non_historical(Kratos.ACCELERATION)])

        # now check
        offset = 0
        for i, node in enumerate(self.model_part.Nodes):
            acceleration = node.GetSolutionStepValue(Kratos.ACCELERATION)
            self.assertEqual(acceleration[0], (data[i + offset, 0] + 2) * 3)
            self.assertEqual(acceleration[1], (data[i + offset, 1] + 2) * 3)
            self.assertEqual(acceleration[2], (data[i + offset, 2] + 2) * 3)

        offset += self.model_part.NumberOfNodes()
        for i, node in enumerate(self.model_part.Nodes):
            acceleration = node.GetValue(Kratos.ACCELERATION)
            self.assertEqual(acceleration[0], (data[i + offset, 0] + 2) * 3)
            self.assertEqual(acceleration[1], (data[i + offset, 1] + 2) * 3)
            self.assertEqual(acceleration[2], (data[i + offset, 2] + 2) * 3)

        offset += self.model_part.NumberOfNodes()
        for i, condition in enumerate(self.model_part.Conditions):
            acceleration = condition.GetValue(Kratos.ACCELERATION)
            self.assertEqual(acceleration[0], (data[i + offset, 0] + 2) * 3)
            self.assertEqual(acceleration[1], (data[i + offset, 1] + 2) * 3)
            self.assertEqual(acceleration[2], (data[i + offset, 2] + 2) * 3)

        offset += self.model_part.NumberOfConditions()
        for i, element in enumerate(self.model_part.Elements):
            acceleration = element.GetValue(Kratos.ACCELERATION)
            self.assertEqual(acceleration[0], (data[i + offset, 0] + 2) * 3)
            self.assertEqual(acceleration[1], (data[i + offset, 1] + 2) * 3)
            self.assertEqual(acceleration[2], (data[i + offset, 2] + 2) * 3)

    def test_Move(self):
        a = Kratos.Expression.NodalExpression(self.model_part)
        b = Kratos.Expression.NodalExpression(self.model_part)
        c = Kratos.Expression.ConditionExpression(self.model_part)
        d = Kratos.Expression.ElementExpression(self.model_part)

        collective = KratosOA.CollectiveExpression([a, b, c, d])

        data = numpy.arange(1, self.model_part.NumberOfNodes() * 4 + self.model_part.NumberOfConditions() * 3 + self.model_part.NumberOfElements() + 1, 1.0)

        # here we move data. The life time is not managed by the collective. If moved data is destroyed, then use of collective can seg fault.
        KratosOA.CollectiveExpressionIO.Move(collective, data, [[3], [], [3], []])

        # do some calculations on all container expressions. these are lazy expressions, hence light wieght operations.
        collective += 2
        collective *= 3

        ceio = KratosOA.CollectiveExpressionIO
        historical = ceio.HistoricalVariable
        non_historical = ceio.NonHistoricalVariable
        KratosOA.CollectiveExpressionIO.Write(collective, [historical(Kratos.ACCELERATION), non_historical(Kratos.DENSITY), non_historical(Kratos.ACCELERATION), non_historical(Kratos.DENSITY)])

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
        KratosOA.CollectiveExpressionIO.Write(collective, [historical(Kratos.ACCELERATION), non_historical(Kratos.DENSITY), non_historical(Kratos.ACCELERATION), non_historical(Kratos.DENSITY)])

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

    def test_Move2(self):
        a = Kratos.Expression.NodalExpression(self.model_part)
        b = Kratos.Expression.NodalExpression(self.model_part)
        c = Kratos.Expression.ConditionExpression(self.model_part)
        d = Kratos.Expression.ElementExpression(self.model_part)

        collective = KratosOA.CollectiveExpression([a, b, c, d])

        data = numpy.arange(1, self.model_part.NumberOfNodes() * 6 + self.model_part.NumberOfConditions() * 3 + self.model_part.NumberOfElements() * 3 + 1, 1.0)
        data = data.reshape((self.model_part.NumberOfNodes() * 2 + self.model_part.NumberOfConditions() + self.model_part.NumberOfElements(), 3))

        # here we move data. The life time is not managed by the collective. If moved data is destroyed, then use of collective can seg fault.
        KratosOA.CollectiveExpressionIO.Move(collective, data)

        # do some calculations on all container expressions. these are lazy expressions, hence light wieght operations.
        collective += 2
        collective *= 3

        ceio = KratosOA.CollectiveExpressionIO
        historical = ceio.HistoricalVariable
        non_historical = ceio.NonHistoricalVariable
        KratosOA.CollectiveExpressionIO.Write(collective, [historical(Kratos.ACCELERATION), non_historical(Kratos.ACCELERATION), non_historical(Kratos.ACCELERATION), non_historical(Kratos.ACCELERATION)])

        # now check for initial values
        offset = 0
        for i, node in enumerate(self.model_part.Nodes):
            acceleration = node.GetSolutionStepValue(Kratos.ACCELERATION)
            self.assertEqual(acceleration[0], (data[i + offset, 0] + 2) * 3)
            self.assertEqual(acceleration[1], (data[i + offset, 1] + 2) * 3)
            self.assertEqual(acceleration[2], (data[i + offset, 2] + 2) * 3)

        offset += self.model_part.NumberOfNodes()
        for i, node in enumerate(self.model_part.Nodes):
            acceleration = node.GetValue(Kratos.ACCELERATION)
            self.assertEqual(acceleration[0], (data[i + offset, 0] + 2) * 3)
            self.assertEqual(acceleration[1], (data[i + offset, 1] + 2) * 3)
            self.assertEqual(acceleration[2], (data[i + offset, 2] + 2) * 3)

        offset += self.model_part.NumberOfNodes()
        for i, condition in enumerate(self.model_part.Conditions):
            acceleration = condition.GetValue(Kratos.ACCELERATION)
            self.assertEqual(acceleration[0], (data[i + offset, 0] + 2) * 3)
            self.assertEqual(acceleration[1], (data[i + offset, 1] + 2) * 3)
            self.assertEqual(acceleration[2], (data[i + offset, 2] + 2) * 3)

        offset += self.model_part.NumberOfConditions()
        for i, element in enumerate(self.model_part.Elements):
            acceleration = element.GetValue(Kratos.ACCELERATION)
            self.assertEqual(acceleration[0], (data[i + offset, 0] + 2) * 3)
            self.assertEqual(acceleration[1], (data[i + offset, 1] + 2) * 3)
            self.assertEqual(acceleration[2], (data[i + offset, 2] + 2) * 3)

        # now change the numpy array to check whether we are still referening to data from numpy arraydata
        data += 1.0
        KratosOA.CollectiveExpressionIO.Write(collective, [historical(Kratos.ACCELERATION), non_historical(Kratos.ACCELERATION), non_historical(Kratos.ACCELERATION), non_historical(Kratos.ACCELERATION)])

        # now check for changed values
        offset = 0
        for i, node in enumerate(self.model_part.Nodes):
            acceleration = node.GetSolutionStepValue(Kratos.ACCELERATION)
            self.assertEqual(acceleration[0], (data[i + offset, 0] + 2) * 3)
            self.assertEqual(acceleration[1], (data[i + offset, 1] + 2) * 3)
            self.assertEqual(acceleration[2], (data[i + offset, 2] + 2) * 3)

        offset += self.model_part.NumberOfNodes()
        for i, node in enumerate(self.model_part.Nodes):
            acceleration = node.GetValue(Kratos.ACCELERATION)
            self.assertEqual(acceleration[0], (data[i + offset, 0] + 2) * 3)
            self.assertEqual(acceleration[1], (data[i + offset, 1] + 2) * 3)
            self.assertEqual(acceleration[2], (data[i + offset, 2] + 2) * 3)

        offset += self.model_part.NumberOfNodes()
        for i, condition in enumerate(self.model_part.Conditions):
            acceleration = condition.GetValue(Kratos.ACCELERATION)
            self.assertEqual(acceleration[0], (data[i + offset, 0] + 2) * 3)
            self.assertEqual(acceleration[1], (data[i + offset, 1] + 2) * 3)
            self.assertEqual(acceleration[2], (data[i + offset, 2] + 2) * 3)

        offset += self.model_part.NumberOfConditions()
        for i, element in enumerate(self.model_part.Elements):
            acceleration = element.GetValue(Kratos.ACCELERATION)
            self.assertEqual(acceleration[0], (data[i + offset, 0] + 2) * 3)
            self.assertEqual(acceleration[1], (data[i + offset, 1] + 2) * 3)
            self.assertEqual(acceleration[2], (data[i + offset, 2] + 2) * 3)

    def test_ReadMoveErrors(self):
        a = Kratos.Expression.NodalExpression(self.model_part)
        b = Kratos.Expression.NodalExpression(self.model_part)
        c = Kratos.Expression.ConditionExpression(self.model_part)
        d = Kratos.Expression.ElementExpression(self.model_part)

        collective = KratosOA.CollectiveExpression([a, b, c, d])

        total_entities = self.model_part.NumberOfNodes() * 4 + self.model_part.NumberOfConditions() * 3 + self.model_part.NumberOfElements() + 1

        with self.assertRaises(TypeError):
            numpy_array = numpy.arange(0, total_entities)
            KratosOA.CollectiveExpressionIO.Read(collective, numpy_array, [[3], [], [3], []])

        with self.assertRaises(TypeError):
            numpy_array = numpy.arange(0, total_entities, dtype=numpy.float32)
            KratosOA.CollectiveExpressionIO.Read(collective, numpy_array, [[3], [], [3], []])

        with self.assertRaises(TypeError):
            numpy_array = numpy.arange(0, total_entities, dtype=numpy.int32)
            KratosOA.CollectiveExpressionIO.Read(collective, numpy_array, [[3], [], [3], []])

        with self.assertRaises(TypeError):
            numpy_array = numpy.arange(0, total_entities, dtype=numpy.int64)
            KratosOA.CollectiveExpressionIO.Read(collective, numpy_array, [[3], [], [3], []])

        with self.assertRaises(TypeError):
            numpy_array = numpy.arange(0, total_entities)
            KratosOA.CollectiveExpressionIO.Move(collective, numpy_array, [[3], [], [3], []])

        with self.assertRaises(TypeError):
            numpy_array = numpy.arange(0, total_entities, dtype=numpy.float32)
            KratosOA.CollectiveExpressionIO.Move(collective, numpy_array, [[3], [], [3], []])

        with self.assertRaises(TypeError):
            numpy_array = numpy.arange(0, total_entities, dtype=numpy.int32)
            KratosOA.CollectiveExpressionIO.Move(collective, numpy_array, [[3], [], [3], []])

        with self.assertRaises(TypeError):
            numpy_array = numpy.arange(0, total_entities, dtype=numpy.int64)
            KratosOA.CollectiveExpressionIO.Move(collective, numpy_array, [[3], [], [3], []])

    def test_CollectiveNegation(self):
        a = Kratos.Expression.NodalExpression(self.model_part)
        b = Kratos.Expression.NodalExpression(self.model_part)
        c = Kratos.Expression.ConditionExpression(self.model_part)
        d = Kratos.Expression.ElementExpression(self.model_part)

        number_of_nodes = self.model_part.NumberOfNodes()
        number_of_conditions = self.model_part.NumberOfConditions()
        number_of_elements = self.model_part.NumberOfElements()
        d_a = numpy.arange(1, number_of_nodes * 3 + 1, dtype=numpy.float64).reshape((number_of_nodes, 3))
        d_b = numpy.arange(1, number_of_nodes * 4 + 1, dtype=numpy.float64).reshape((number_of_nodes, 2, 2))
        d_c = numpy.arange(1, number_of_conditions * 6 + 1, dtype=numpy.float64).reshape((number_of_conditions, 2, 3))
        d_d = numpy.arange(1, number_of_elements * 9 + 1, dtype=numpy.float64).reshape((number_of_elements, 3, 3))

        Kratos.Expression.CArrayExpressionIO.Read(a, d_a)
        Kratos.Expression.CArrayExpressionIO.Read(b, d_b)
        Kratos.Expression.CArrayExpressionIO.Read(c, d_c)
        Kratos.Expression.CArrayExpressionIO.Read(d, d_d)

        c_a = KratosOA.CollectiveExpression([a, b, c, d])
        c_b = KratosOA.CollectiveExpression([-a, b * 2, c + 3, d / 2])

        d = -c_a
        self.assertVectorAlmostEqual(c_a.Evaluate(), numpy.concatenate([d_a.flatten(), d_b.flatten(), d_c.flatten(), d_d.flatten()]).flatten())
        self.assertVectorAlmostEqual(c_b.Evaluate(), numpy.concatenate([-d_a.flatten(), d_b.flatten() * 2, d_c.flatten() + 3, d_d.flatten() / 2]).flatten())
        self.assertVectorAlmostEqual(d.Evaluate(), -numpy.concatenate([d_a.flatten(), d_b.flatten(), d_c.flatten(), d_d.flatten()]).flatten())

        d -= c_b
        self.assertVectorAlmostEqual(d.Evaluate(), -numpy.concatenate([d_a.flatten(), d_b.flatten(), d_c.flatten(), d_d.flatten()]).flatten() - numpy.concatenate([-d_a.flatten(), d_b.flatten() * 2, d_c.flatten() + 3, d_d.flatten() / 2]).flatten())
        self.assertVectorAlmostEqual(c_a.Evaluate(), numpy.concatenate([d_a.flatten(), d_b.flatten(), d_c.flatten(), d_d.flatten()]).flatten())
        self.assertVectorAlmostEqual(c_b.Evaluate(), numpy.concatenate([-d_a.flatten(), d_b.flatten() * 2, d_c.flatten() + 3, d_d.flatten() / 2]).flatten())

    def test_CollectiveAddition(self):
        a = Kratos.Expression.NodalExpression(self.model_part)
        b = Kratos.Expression.NodalExpression(self.model_part)
        c = Kratos.Expression.ConditionExpression(self.model_part)
        d = Kratos.Expression.ElementExpression(self.model_part)

        number_of_nodes = self.model_part.NumberOfNodes()
        number_of_conditions = self.model_part.NumberOfConditions()
        number_of_elements = self.model_part.NumberOfElements()
        d_a = numpy.arange(1, number_of_nodes * 3 + 1, dtype=numpy.float64).reshape((number_of_nodes, 3))
        d_b = numpy.arange(1, number_of_nodes * 4 + 1, dtype=numpy.float64).reshape((number_of_nodes, 2, 2))
        d_c = numpy.arange(1, number_of_conditions * 6 + 1, dtype=numpy.float64).reshape((number_of_conditions, 2, 3))
        d_d = numpy.arange(1, number_of_elements * 9 + 1, dtype=numpy.float64).reshape((number_of_elements, 3, 3))

        Kratos.Expression.CArrayExpressionIO.Read(a, d_a)
        Kratos.Expression.CArrayExpressionIO.Read(b, d_b)
        Kratos.Expression.CArrayExpressionIO.Read(c, d_c)
        Kratos.Expression.CArrayExpressionIO.Read(d, d_d)

        c_a = KratosOA.CollectiveExpression([a, b, c, d])
        c_b = KratosOA.CollectiveExpression([-a, b * 2, c + 3, d / 2])

        d = c_a.Clone()
        d += c_b * 3
        self.assertVectorAlmostEqual(d.Evaluate(), numpy.concatenate([d_a.flatten(), d_b.flatten(), d_c.flatten(), d_d.flatten()]).flatten() + numpy.concatenate([-d_a.flatten(), d_b.flatten() * 2, d_c.flatten() + 3, d_d.flatten() / 2]).flatten() * 3)
        self.assertVectorAlmostEqual(c_a.Evaluate(), numpy.concatenate([d_a.flatten(), d_b.flatten(), d_c.flatten(), d_d.flatten()]).flatten())
        self.assertVectorAlmostEqual(c_b.Evaluate(), numpy.concatenate([-d_a.flatten(), d_b.flatten() * 2, d_c.flatten() + 3, d_d.flatten() / 2]).flatten())

    def test_CollectiveMultiplication(self):
        a = Kratos.Expression.NodalExpression(self.model_part)
        b = Kratos.Expression.NodalExpression(self.model_part)
        c = Kratos.Expression.ConditionExpression(self.model_part)
        d = Kratos.Expression.ElementExpression(self.model_part)

        number_of_nodes = self.model_part.NumberOfNodes()
        number_of_conditions = self.model_part.NumberOfConditions()
        number_of_elements = self.model_part.NumberOfElements()
        d_a = numpy.arange(1, number_of_nodes * 3 + 1, dtype=numpy.float64).reshape((number_of_nodes, 3))
        d_b = numpy.arange(1, number_of_nodes * 4 + 1, dtype=numpy.float64).reshape((number_of_nodes, 2, 2))
        d_c = numpy.arange(1, number_of_conditions * 6 + 1, dtype=numpy.float64).reshape((number_of_conditions, 2, 3))
        d_d = numpy.arange(1, number_of_elements * 9 + 1, dtype=numpy.float64).reshape((number_of_elements, 3, 3))

        Kratos.Expression.CArrayExpressionIO.Read(a, d_a)
        Kratos.Expression.CArrayExpressionIO.Read(b, d_b)
        Kratos.Expression.CArrayExpressionIO.Read(c, d_c)
        Kratos.Expression.CArrayExpressionIO.Read(d, d_d)

        c_a = KratosOA.CollectiveExpression([a, b, c, d])
        c_b = KratosOA.CollectiveExpression([-a, b * 2, c + 3, d / 2])

        d = c_a.Clone()
        d *= c_b
        self.assertVectorAlmostEqual(d.Evaluate(), numpy.concatenate([d_a.flatten(), d_b.flatten(), d_c.flatten(), d_d.flatten()]).flatten() * numpy.concatenate([-d_a.flatten(), d_b.flatten() * 2, d_c.flatten() + 3, d_d.flatten() / 2]).flatten())
        self.assertVectorAlmostEqual(c_a.Evaluate(), numpy.concatenate([d_a.flatten(), d_b.flatten(), d_c.flatten(), d_d.flatten()]).flatten())
        self.assertVectorAlmostEqual(c_b.Evaluate(), numpy.concatenate([-d_a.flatten(), d_b.flatten() * 2, d_c.flatten() + 3, d_d.flatten() / 2]).flatten())

    def test_CollectiveDivision(self):
        a = Kratos.Expression.NodalExpression(self.model_part)
        b = Kratos.Expression.NodalExpression(self.model_part)
        c = Kratos.Expression.ConditionExpression(self.model_part)
        d = Kratos.Expression.ElementExpression(self.model_part)

        number_of_nodes = self.model_part.NumberOfNodes()
        number_of_conditions = self.model_part.NumberOfConditions()
        number_of_elements = self.model_part.NumberOfElements()
        d_a = numpy.arange(1, number_of_nodes * 3 + 1, dtype=numpy.float64).reshape((number_of_nodes, 3))
        d_b = numpy.arange(1, number_of_nodes * 4 + 1, dtype=numpy.float64).reshape((number_of_nodes, 2, 2))
        d_c = numpy.arange(1, number_of_conditions * 6 + 1, dtype=numpy.float64).reshape((number_of_conditions, 2, 3))
        d_d = numpy.arange(1, number_of_elements * 9 + 1, dtype=numpy.float64).reshape((number_of_elements, 3, 3))

        Kratos.Expression.CArrayExpressionIO.Read(a, d_a)
        Kratos.Expression.CArrayExpressionIO.Read(b, d_b)
        Kratos.Expression.CArrayExpressionIO.Read(c, d_c)
        Kratos.Expression.CArrayExpressionIO.Read(d, d_d)

        c_a = KratosOA.CollectiveExpression([a, b, c, d])
        c_b = KratosOA.CollectiveExpression([-a, b * 2, c + 3, d / 2])

        d = c_a.Clone()
        d /= c_b
        self.assertVectorAlmostEqual(d.Evaluate(), numpy.concatenate([d_a.flatten(), d_b.flatten(), d_c.flatten(), d_d.flatten()]).flatten() / numpy.concatenate([-d_a.flatten(), d_b.flatten() * 2, d_c.flatten() + 3, d_d.flatten() / 2]).flatten())
        self.assertVectorAlmostEqual(c_a.Evaluate(), numpy.concatenate([d_a.flatten(), d_b.flatten(), d_c.flatten(), d_d.flatten()]).flatten())
        self.assertVectorAlmostEqual(c_b.Evaluate(), numpy.concatenate([-d_a.flatten(), d_b.flatten() * 2, d_c.flatten() + 3, d_d.flatten() / 2]).flatten())


if __name__ == "__main__":
    Kratos.Tester.SetVerbosity(Kratos.Tester.Verbosity.PROGRESS)  # TESTS_OUTPUTS
    kratos_unittest.main()