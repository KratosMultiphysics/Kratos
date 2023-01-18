
import KratosMultiphysics as Kratos

# Import KratosUnittest
import KratosMultiphysics.KratosUnittest as kratos_unittest

from KratosMultiphysics.OptimizationApplication.utilities.container_data import ContainerData

class TestContainerData(kratos_unittest.TestCase):
    def setUp(self):
        self.model = Kratos.Model()
        self.model_part = self.model.CreateModelPart("test")
        self.model_part.AddNodalSolutionStepVariable(Kratos.DENSITY)
        self.model_part.AddNodalSolutionStepVariable(Kratos.ACCELERATION)

        number_of_nodes = 10
        for id in range(1, number_of_nodes + 1):
            node = self.model_part.CreateNewNode(id, id, id+1, id+2)
            node.SetValue(Kratos.PRESSURE, id+3)
            node.SetValue(Kratos.VELOCITY, Kratos.Array3([id+3, id+4, id+5]))

        number_of_conditions = 11
        for id in range(1, number_of_conditions + 1):
            properties = self.model_part.CreateNewProperties(id)
            properties.SetValue(Kratos.PRESSURE, id+400)
            properties.SetValue(Kratos.VELOCITY, Kratos.Array3([id+500, id+600, id+700]))
            condition = self.model_part.CreateNewCondition("LineCondition2D2N", id + 1, [(id % number_of_nodes) + 1, ((id + 1) % number_of_nodes) + 1 ], properties)
            condition.SetValue(Kratos.PRESSURE, id+4)
            condition.SetValue(Kratos.VELOCITY, Kratos.Array3([id+5, id+6, id+7]))

        number_of_elements = 12
        for id in range(1, number_of_elements + 1):
            properties = self.model_part.CreateNewProperties(id + number_of_conditions)
            properties.SetValue(Kratos.PRESSURE, id+500)
            properties.SetValue(Kratos.VELOCITY, Kratos.Array3([id+600, id+700, id+800]))
            element = self.model_part.CreateNewElement("Element2D3N", id + 2, [(id % number_of_nodes) + 1, ((id + 1) % number_of_nodes) + 1, ((id + 2) % number_of_nodes) + 1 ], properties)
            element.SetValue(Kratos.PRESSURE, id+5)
            element.SetValue(Kratos.VELOCITY, Kratos.Array3([id+6, id+7, id+8]))


    def test_ContainerDataAdd(self):
        a = ContainerData(self.model_part, ContainerData.ContainerEnum.ELEMENT_PROPERTIES)
        b = ContainerData(self.model_part, ContainerData.ContainerEnum.ELEMENT_PROPERTIES)

        a.ReadDataFromContainer(Kratos.PRESSURE)
        b.ReadDataFromContainer(Kratos.PRESSURE)
        c = a + b

        c.AssignDataToContainer(Kratos.DENSITY)
        for element in self.model_part.Elements:
            self.assertEqual(element.Properties[Kratos.DENSITY], element.Properties[Kratos.PRESSURE] * 2)

    def test_ContainerDataSubstract(self):
        a = ContainerData(self.model_part, ContainerData.ContainerEnum.ELEMENT_PROPERTIES)
        b = ContainerData(self.model_part, ContainerData.ContainerEnum.ELEMENT_PROPERTIES)

        a.ReadDataFromContainer(Kratos.PRESSURE)
        b.ReadDataFromContainer(Kratos.PRESSURE)
        c = a * 4 - b

        c.AssignDataToContainer(Kratos.DENSITY)
        for element in self.model_part.Elements:
            self.assertEqual(element.Properties[Kratos.DENSITY], element.Properties[Kratos.PRESSURE] * 3)

    def test_ContainerDataMultiply(self):
        a = ContainerData(self.model_part, ContainerData.ContainerEnum.ELEMENT_PROPERTIES)

        a.ReadDataFromContainer(Kratos.PRESSURE)
        c = a * 4

        c.AssignDataToContainer(Kratos.DENSITY)
        for element in self.model_part.Elements:
            self.assertEqual(element.Properties[Kratos.DENSITY], element.Properties[Kratos.PRESSURE] * 4)

    def test_ContainerDataDivide(self):
        a = ContainerData(self.model_part, ContainerData.ContainerEnum.ELEMENT_PROPERTIES)

        a.ReadDataFromContainer(Kratos.PRESSURE)
        c = a / 2

        c.AssignDataToContainer(Kratos.DENSITY)
        for element in self.model_part.Elements:
            self.assertEqual(element.Properties[Kratos.DENSITY], element.Properties[Kratos.PRESSURE] / 2)

    def test_IsSameContainer(self):
        a = ContainerData(self.model_part, ContainerData.ContainerEnum.ELEMENT_PROPERTIES)

        a.ReadDataFromContainer(Kratos.PRESSURE)
        b = a.Clone()
        b = b * 2
        b.AssignDataToContainer(Kratos.DENSITY)

        self.assertTrue(a.IsSameContainer(b))
        self.assertVectorAlmostEqual(a.GetData(), b.GetData() / 2, 9)

if __name__ == "__main__":
    Kratos.Tester.SetVerbosity(Kratos.Tester.Verbosity.PROGRESS)  # TESTS_OUTPUTS
    kratos_unittest.main()