
import KratosMultiphysics as Kratos
import KratosMultiphysics.OptimizationApplication as KratosOA

# Import KratosUnittest
import KratosMultiphysics.KratosUnittest as kratos_unittest

class TestOptimizationVariableUtils(kratos_unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.model = Kratos.Model()
        cls.model_part = cls.model.CreateModelPart("test")

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

        cls.utils = KratosOA.OptimizationVariableUtils

    def test_GetContainerIds(self):
        # check for nodes
        self.assertEqual(self.utils.GetContainerIds(self.model_part.Nodes), [node.Id for node in self.model_part.Nodes])

        # check for conditions
        self.assertEqual(self.utils.GetContainerIds(self.model_part.Conditions), [condition.Id for condition in self.model_part.Conditions])

        # check for elements
        self.assertEqual(self.utils.GetContainerIds(self.model_part.Elements), [element.Id for element in self.model_part.Elements])

    def test_GetContainerVariableToMatrix(self):
        values = Kratos.Matrix()
        self.utils.GetContainerVariableToMatrix(self.model_part.Nodes, Kratos.VELOCITY, values)
        self.__CheckMatrix(values, self.utils.GetContainerIds(self.model_part.Nodes), lambda m, id: m.GetNode(id).GetValue(Kratos.VELOCITY))

        self.utils.GetContainerVariableToMatrix(self.model_part.Conditions, Kratos.VELOCITY, values)
        self.__CheckMatrix(values, self.utils.GetContainerIds(self.model_part.Conditions), lambda m, id: m.GetCondition(id).GetValue(Kratos.VELOCITY))

        self.utils.GetContainerVariableToMatrix(self.model_part.Elements, Kratos.VELOCITY, values)
        self.__CheckMatrix(values, self.utils.GetContainerIds(self.model_part.Elements), lambda m, id: m.GetElement(id).GetValue(Kratos.VELOCITY))

    def test_GetContainerVariableToVector(self):
        values = Kratos.Vector()
        self.utils.GetContainerVariableToVector(self.model_part.Nodes, Kratos.PRESSURE, values)
        self.__CheckVector(values, self.utils.GetContainerIds(self.model_part.Nodes), lambda m, id: m.GetNode(id).GetValue(Kratos.PRESSURE))

        self.utils.GetContainerVariableToVector(self.model_part.Conditions, Kratos.PRESSURE, values)
        self.__CheckVector(values, self.utils.GetContainerIds(self.model_part.Conditions), lambda m, id: m.GetCondition(id).GetValue(Kratos.PRESSURE))

        self.utils.GetContainerVariableToVector(self.model_part.Elements, Kratos.PRESSURE, values)
        self.__CheckVector(values, self.utils.GetContainerIds(self.model_part.Elements), lambda m, id: m.GetElement(id).GetValue(Kratos.PRESSURE))

    def test_GetContainerPropertiesVariableToVector(self):
        values = Kratos.Vector()
        self.utils.GetContainerPropertiesVariableToVector(self.model_part.Conditions, Kratos.PRESSURE, values)
        self.__CheckVector(values, self.utils.GetContainerIds(self.model_part.Conditions), lambda m, id: m.GetCondition(id).Properties[Kratos.PRESSURE])

        self.utils.GetContainerPropertiesVariableToVector(self.model_part.Elements, Kratos.PRESSURE, values)
        self.__CheckVector(values, self.utils.GetContainerIds(self.model_part.Elements), lambda m, id: m.GetElement(id).Properties[Kratos.PRESSURE])

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

    def __CheckMatrix(self, matrix: Kratos.Matrix, ids, get_entity_value_method):
        self.assertEqual(matrix.Size1(), len(ids))
        for i, id in enumerate(ids):
            v = get_entity_value_method(self.model_part, id)
            for j in range(matrix.Size2()):
                self.assertEqual(matrix[i, j], v[j])

    def __CheckVector(self, vector: Kratos.Vector, ids, get_entity_value_method):
        self.assertEqual(vector.Size(), len(ids))
        for i, id in enumerate(ids):
            v = get_entity_value_method(self.model_part, id)
            self.assertEqual(v, vector[i])

if __name__ == "__main__":
    Kratos.Tester.SetVerbosity(Kratos.Tester.Verbosity.PROGRESS)  # TESTS_OUTPUTS
    kratos_unittest.main()