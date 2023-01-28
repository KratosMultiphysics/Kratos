
import KratosMultiphysics as Kratos
import KratosMultiphysics.OptimizationApplication as KratosOA

# Import KratosUnittest
import KratosMultiphysics.KratosUnittest as kratos_unittest

class TestContainerData(kratos_unittest.TestCase):
    def setUp(self):
        self.model = Kratos.Model()
        self.model_part = self.model.CreateModelPart("test")
        self.model_part.AddNodalSolutionStepVariable(Kratos.DENSITY)
        self.model_part.AddNodalSolutionStepVariable(Kratos.ACCELERATION)
        self.model_part.ProcessInfo[Kratos.DOMAIN_SIZE] = 3

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
        a = KratosOA.ContainerData(self.model_part, KratosOA.ContainerDataType.ElementProperties)
        b = KratosOA.ContainerData(self.model_part, KratosOA.ContainerDataType.ElementProperties)

        a.ReadDataFromContainerVariable(Kratos.VELOCITY)
        b.ReadDataFromContainerVariable(Kratos.VELOCITY)

        c = a + b

        c.AssignDataToContainerVariable(Kratos.ACCELERATION)
        for element in c.GetContainer():
            self.assertVectorAlmostEqual(element.Properties[Kratos.ACCELERATION], element.Properties[Kratos.VELOCITY] * 2, 12)

    def test_ContainerDataSubstract(self):
        a = KratosOA.ContainerData(self.model_part, KratosOA.ContainerDataType.ElementProperties)
        b = KratosOA.ContainerData(self.model_part, KratosOA.ContainerDataType.ElementProperties)

        a.ReadDataFromContainerVariable(Kratos.VELOCITY)
        b.ReadDataFromContainerVariable(Kratos.VELOCITY)
        c = a * 4 - b

        c.AssignDataToContainerVariable(Kratos.ACCELERATION)
        for element in c.GetContainer():
            self.assertVectorAlmostEqual(element.Properties[Kratos.ACCELERATION], element.Properties[Kratos.VELOCITY] * 3, 12)

    def test_ContainerDataMultiply(self):
        a = KratosOA.ContainerData(self.model_part, KratosOA.ContainerDataType.ElementProperties)

        a.ReadDataFromContainerVariable(Kratos.VELOCITY)
        c = a * 4

        c.AssignDataToContainerVariable(Kratos.ACCELERATION)
        for element in c.GetContainer():
            self.assertVectorAlmostEqual(element.Properties[Kratos.ACCELERATION], element.Properties[Kratos.VELOCITY] * 4, 12)

    def test_ContainerDataDivide(self):
        a = KratosOA.ContainerData(self.model_part, KratosOA.ContainerDataType.ElementProperties)

        a.ReadDataFromContainerVariable(Kratos.VELOCITY)
        c = a / 2

        c.AssignDataToContainerVariable(Kratos.ACCELERATION)
        for element in c.GetContainer():
            self.assertVectorAlmostEqual(element.Properties[Kratos.ACCELERATION], element.Properties[Kratos.VELOCITY] / 2, 12)

    def test_IsSameContainer(self):
        a = KratosOA.ContainerData(self.model_part, KratosOA.ContainerDataType.ElementProperties)

        a.ReadDataFromContainerVariable(Kratos.VELOCITY)
        b = a.Clone()
        b = b * 2
        b.AssignDataToContainerVariable(Kratos.ACCELERATION)

        self.assertTrue(a.IsSameContainer(b))
        self.assertVectorAlmostEqual(a.GetData(), b.GetData() / 2, 9)

    def test_SetDataForContainerVariable(self):
        a = KratosOA.ContainerData(self.model_part, KratosOA.ContainerDataType.ElementProperties)
        a.SetDataForContainerVariable(Kratos.VELOCITY, Kratos.Array3([1, 2, 3]))
        a.AssignDataToContainerVariable(Kratos.ACCELERATION)
        for element in a.GetContainer():
            self.assertVectorAlmostEqual(element.Properties[Kratos.ACCELERATION], Kratos.Array3([1, 2, 3]), 12)

    def test_Clone(self):
        a = KratosOA.ContainerData(self.model_part, KratosOA.ContainerDataType.ElementProperties)

        a.ReadDataFromContainerVariable(Kratos.VELOCITY)

        b = a.Clone()
        b.SetDataForContainerVariable(Kratos.VELOCITY, Kratos.Array3([10, 11, 12]))

        a.AssignDataToContainerVariable(Kratos.ACCELERATION)
        for element in self.model_part.Elements:
            self.assertVectorAlmostEqual(element.Properties[Kratos.ACCELERATION], element.Properties[Kratos.VELOCITY], 12)

        b.AssignDataToContainerVariable(Kratos.ACCELERATION)
        for element in a.GetContainer():
            self.assertVectorAlmostEqual(element.Properties[Kratos.ACCELERATION], Kratos.Array3([10, 11, 12]), 12)

if __name__ == "__main__":
    Kratos.Tester.SetVerbosity(Kratos.Tester.Verbosity.PROGRESS)  # TESTS_OUTPUTS
    kratos_unittest.main()