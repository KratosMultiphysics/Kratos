
import KratosMultiphysics as Kratos
import KratosMultiphysics.OptimizationApplication as KratosOA

# Import KratosUnittest
import KratosMultiphysics.KratosUnittest as kratos_unittest



class TestEntitySpecificPropertiesProcess(kratos_unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.model = Kratos.Model()
        cls.model_part = cls.model.CreateModelPart("test")
        properties = cls.model_part.CreateNewProperties(1)

        properties[Kratos.DENSITY] = 10.0
        properties[Kratos.DIAMETER] = 15.0
        properties[Kratos.PRESSURE] = -5.0
        properties[Kratos.VELOCITY] = Kratos.Array3([1, 2, 3])

        number_of_entities = 10
        for i in range(1, number_of_entities + 1):
            cls.model_part.CreateNewNode(i, 0, 0, 0)

        for i in range(1, number_of_entities + 1):
            cls.model_part.CreateNewCondition("LineCondition2D2N", i, [(i % number_of_entities) + 1, ((i + 1) % number_of_entities) + 1], properties)

        for i in range(1, number_of_entities + 1):
            cls.model_part.CreateNewElement("Element2D3N", i, [(i % number_of_entities) + 1, ((i + 1) % number_of_entities) + 1, ((i + 2) % number_of_entities) + 1], properties)


    def test_EntitySpecificPropertiesProcess(self):
        parameters = Kratos.Parameters("""{
            "model_part_name": "test",
            "variables_list" : ["DENSITY", "PRESSURE"],
            "echo_level"     : 0
        }""")

        process = KratosOA.EntitySpecificPropertiesProcess(self.model, parameters)

        # before the process, following variables shoudn't be in the element/condition data value container
        for element in self.model_part.Elements:
            self.assertFalse(element.Has(Kratos.DENSITY))
            self.assertFalse(element.Has(Kratos.PRESSURE))

        for condition in self.model_part.Conditions:
            self.assertFalse(condition.Has(Kratos.DENSITY))
            self.assertFalse(condition.Has(Kratos.PRESSURE))

        process.ExecuteInitializeSolutionStep()

        # after the first iteration, following variables should be in the element/condition data value container with the
        # values from the properties
        for element in self.model_part.Elements:
            self.assertEqual(element.GetValue(Kratos.DENSITY), 10.0)
            self.assertEqual(element.GetValue(Kratos.PRESSURE), -5.0)

        for condition in self.model_part.Conditions:
            self.assertEqual(condition.GetValue(Kratos.DENSITY), 10.0)
            self.assertEqual(condition.GetValue(Kratos.PRESSURE), -5.0)

        # update the element/condition data container values
        for element in self.model_part.Elements:
            element.SetValue(Kratos.DENSITY, element.GetValue(Kratos.DENSITY) + 10.0)
            element.SetValue(Kratos.PRESSURE, element.GetValue(Kratos.PRESSURE) - 5.0)

        for condition in self.model_part.Conditions:
            condition.SetValue(Kratos.DENSITY, condition.GetValue(Kratos.DENSITY) + 30.0)
            condition.SetValue(Kratos.PRESSURE, condition.GetValue(Kratos.PRESSURE) - 25.0)

        process.ExecuteInitializeSolutionStep()

        # now checking for updated properties
        for element in self.model_part.Elements:
            self.assertEqual(element.Properties[Kratos.DENSITY], 20.0)
            self.assertEqual(element.Properties[Kratos.PRESSURE], -10.0)
            self.assertEqual(element.Properties[Kratos.DIAMETER], 15.0)
            self.assertVectorAlmostEqual(element.Properties[Kratos.VELOCITY], Kratos.Array3([1, 2, 3]), 9)

        for condition in self.model_part.Conditions:
            self.assertEqual(condition.Properties[Kratos.DENSITY], 40.0)
            self.assertEqual(condition.Properties[Kratos.PRESSURE], -30.0)
            self.assertEqual(condition.Properties[Kratos.DIAMETER], 15.0)
            self.assertVectorAlmostEqual(condition.Properties[Kratos.VELOCITY], Kratos.Array3([1, 2, 3]), 9)


if __name__ == "__main__":
    Kratos.Tester.SetVerbosity(Kratos.Tester.Verbosity.PROGRESS)  # TESTS_OUTPUTS
    kratos_unittest.main()