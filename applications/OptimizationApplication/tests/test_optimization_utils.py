
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

    def test_ResetModelPartNodalSolutionStepData(self):

        model = Kratos.Model()
        model_part = model.CreateModelPart("test")

        model_part.AddNodalSolutionStepVariable(Kratos.PRESSURE)
        model_part.AddNodalSolutionStepVariable(Kratos.VELOCITY)

        for i in range(50):
            node: Kratos.Node = model_part.CreateNewNode(i + 1, 0.0, 0.0, 0.0)
            node.SetSolutionStepValue(Kratos.PRESSURE, i + 100)
            node.SetSolutionStepValue(Kratos.VELOCITY, Kratos.Array3([i + 200, i + 300, i + 400]))

        for node in model_part.Nodes:
            self.assertEqual(node.GetSolutionStepValue(Kratos.PRESSURE), node.Id + 99)
            self.assertVectorAlmostEqual(node.GetSolutionStepValue(Kratos.VELOCITY), Kratos.Array3([node.Id+199, node.Id+299, node.Id+399]))

        KratosOA.OptimizationUtils.ResetModelPartNodalSolutionStepData(model_part)

        for node in model_part.Nodes:
            self.assertEqual(node.GetSolutionStepValue(Kratos.PRESSURE), 0.0)
            self.assertVectorAlmostEqual(node.GetSolutionStepValue(Kratos.VELOCITY), Kratos.Array3([0.0, 0.0, 0.0]))

if __name__ == "__main__":
    kratos_unittest.main()