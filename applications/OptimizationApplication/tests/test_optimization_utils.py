
import numpy
import KratosMultiphysics as Kratos
import KratosMultiphysics.OptimizationApplication as KratosOA

# Import KratosUnittest
import KratosMultiphysics.KratosUnittest as kratos_unittest

class TestOptimizationUtils(kratos_unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.model = Kratos.Model()
        cls.model_part = cls.model.CreateModelPart("test")
        cls.model_part.CreateSubModelPart("sub_1")
        cls.model_part.CreateSubModelPart("sub_2")
        cls.model_part.CreateSubModelPart("sub_3")
        cls.model_part.AddNodalSolutionStepVariable(Kratos.DENSITY)
        cls.model_part.AddNodalSolutionStepVariable(Kratos.ACCELERATION)

        number_of_nodes = 10
        for id in range(1, number_of_nodes + 1):
            node = cls.model_part.CreateNewNode(id, id, id+1, id+2)
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

        cls.utils = KratosOA.OptimizationUtils

    def test_GetVariableDimension(self):
        self.assertEqual(self.utils.GetVariableDimension(Kratos.PRESSURE, 1), 1)
        self.assertEqual(self.utils.GetVariableDimension(Kratos.VELOCITY, 2), 2)

    def test_IsVariableExistsInAllContainerProperties(self):
        self.assertTrue(self.utils.IsVariableExistsInAllContainerProperties(self.model_part.Conditions, Kratos.PRESSURE, self.model_part.GetCommunicator().GetDataCommunicator()))
        for condition in self.model_part.Conditions:
            if condition.Id != 4:
                condition.Properties[Kratos.DISTANCE] = 1.0
        self.assertFalse(self.utils.IsVariableExistsInAllContainerProperties(self.model_part.Conditions, Kratos.DISTANCE, self.model_part.GetCommunicator().GetDataCommunicator()))

    def test_IsVariableExistsInAtLeastOneContainerProperties(self):
        self.assertFalse(self.utils.IsVariableExistsInAtLeastOneContainerProperties(self.model_part.Conditions, Kratos.DENSITY, self.model_part.GetCommunicator().GetDataCommunicator()))
        self.model_part.GetCondition(4).Properties[Kratos.DENSITY] = 1.0
        self.assertTrue(self.utils.IsVariableExistsInAtLeastOneContainerProperties(self.model_part.Conditions, Kratos.DENSITY, self.model_part.GetCommunicator().GetDataCommunicator()))

    def test_GetComponentWiseModelParts(self):
        params = Kratos.Parameters("""{
            "test.sub_1": [true, false],
            "test.sub_2": [true, true]
        }""")
        values = KratosOA.OptimizationUtils.GetComponentWiseModelParts(self.model, params)
        self.assertEqual(values[0], [self.model.GetModelPart("test.sub_1"), self.model.GetModelPart("test.sub_2")])
        self.assertEqual(values[1], [self.model.GetModelPart("test.sub_2")])

    def test_CreateEntitySpecificPropertiesForContainer(self):
        model = Kratos.Model()
        model_part = model.CreateModelPart("test_model_part")

        def create_properties(properties_id: int) -> Kratos.Properties:
            properties: Kratos.Properties = model_part.CreateNewProperties(properties_id)
            properties.SetValue(Kratos.PRESSURE, properties_id)
            return properties

        properties = create_properties(1)
        properties.AddSubProperties(create_properties(2))
        properties.AddSubProperties(create_properties(3))
        properties.AddSubProperties(create_properties(4))

        properties.GetSubProperties(2).AddSubProperties(create_properties(5))
        properties.GetSubProperties(3).AddSubProperties(create_properties(6))
        properties.GetSubProperties(4).AddSubProperties(create_properties(7))

        properties.GetSubProperties(4).GetSubProperties(7).AddSubProperties(create_properties(8))

        for i in range(50):
            model_part.CreateNewNode(i + 1, 0, 0, 0)
            element: Kratos.Element = model_part.CreateNewElement("Element3D1N", i + 1, [i+1], properties)

        KratosOA.OptimizationUtils.CreateEntitySpecificPropertiesForContainer(model_part, model_part.Elements, False)

        values = numpy.arange(50, dtype=numpy.float64)
        expression = Kratos.Expression.ElementExpression(model_part)
        Kratos.Expression.CArrayExpressionIO.Read(expression, values)
        KratosOA.PropertiesVariableExpressionIO.Write(expression, Kratos.YOUNG_MODULUS)

        for element in model_part.Elements:
            self.assertEqual(element.Properties.Id, element.Id + 8)
            self.assertEqual(element.Properties[Kratos.PRESSURE], 1)
            self.assertEqual(element.Properties[Kratos.YOUNG_MODULUS], element.Id - 1)

    def test_CreateEntitySpecificPropertiesForContainerRecursive(self):
        model = Kratos.Model()
        model_part = model.CreateModelPart("test_model_part")

        def create_properties(properties_id: int) -> Kratos.Properties:
            properties: Kratos.Properties = model_part.CreateNewProperties(properties_id)
            properties.SetValue(Kratos.PRESSURE, properties_id)
            properties.SetValue(Kratos.YOUNG_MODULUS, 10)
            return properties

        properties = create_properties(1)
        properties.AddSubProperties(create_properties(2))
        properties.AddSubProperties(create_properties(3))
        properties.AddSubProperties(create_properties(4))

        properties.GetSubProperties(2).AddSubProperties(create_properties(5))
        properties.GetSubProperties(3).AddSubProperties(create_properties(6))
        properties.GetSubProperties(4).AddSubProperties(create_properties(7))

        properties.GetSubProperties(4).GetSubProperties(7).AddSubProperties(create_properties(8))

        for i in range(50):
            model_part.CreateNewNode(i + 1, 0, 0, 0)
            element: Kratos.Element = model_part.CreateNewElement("Element3D1N", i + 1, [i+1], properties)

        KratosOA.OptimizationUtils.CreateEntitySpecificPropertiesForContainer(model_part, model_part.Elements, True)

        values = numpy.arange(50, dtype=numpy.float64)
        expression = Kratos.Expression.ElementExpression(model_part)
        Kratos.Expression.CArrayExpressionIO.Read(expression, values)
        KratosOA.PropertiesVariableExpressionIO.Write(expression, Kratos.YOUNG_MODULUS)

        def check_properties_value_recursively(properties: Kratos.Properties, value: float):
            self.assertEqual(properties[Kratos.YOUNG_MODULUS], value)
            for sub_properties in properties.GetSubProperties():
                check_properties_value_recursively(sub_properties, value)

        check_id = 9
        for element in model_part.Elements:
            properties = element.Properties

            self.assertEqual(properties.Id, check_id)
            self.assertEqual(properties[Kratos.YOUNG_MODULUS], element.Id - 1)

            sub_properties_itr = iter(properties.GetSubProperties())

            self.assertEqual(next(sub_properties_itr).Id, check_id + 1)
            self.assertEqual(next(iter(properties.GetSubProperties()[check_id + 1].GetSubProperties())).Id, check_id + 2)
            check_properties_value_recursively(properties.GetSubProperties()[check_id + 1], 10)

            self.assertEqual(next(sub_properties_itr).Id, check_id + 3)
            self.assertEqual(next(iter(properties.GetSubProperties()[check_id + 3].GetSubProperties())).Id, check_id + 4)
            check_properties_value_recursively(properties.GetSubProperties()[check_id + 3], 10)

            self.assertEqual(next(sub_properties_itr).Id, check_id + 5)
            self.assertEqual(next(iter(properties.GetSubProperties()[check_id + 5].GetSubProperties())).Id, check_id + 6)
            self.assertEqual(next(iter(properties.GetSubProperties()[check_id + 5].GetSubProperties()[check_id + 6].GetSubProperties())).Id, check_id + 7)
            check_properties_value_recursively(properties.GetSubProperties()[check_id + 5], 10)

            check_id += 8

        KratosOA.OptimizationUtils.UpdatePropertiesVariableWithRootValueRecursively(model_part.Elements, Kratos.YOUNG_MODULUS)
        for element in model_part.Elements:
            check_properties_value_recursively(element.Properties, element.Id - 1)

if __name__ == "__main__":
    kratos_unittest.main()