import math
import KratosMultiphysics as Kratos
import KratosMultiphysics.OptimizationApplication as KratosOA

# Import KratosUnittest
import KratosMultiphysics.KratosUnittest as kratos_unittest

class TestCollectiveVariableDataHolderBase(kratos_unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.model = Kratos.Model()
        cls.model_part = cls.model.CreateModelPart("test")
        cls.model_part.AddNodalSolutionStepVariable(Kratos.DENSITY)
        cls.model_part.AddNodalSolutionStepVariable(Kratos.PRESSURE)
        cls.model_part.AddNodalSolutionStepVariable(Kratos.ACCELERATION)
        cls.model_part.AddNodalSolutionStepVariable(Kratos.VELOCITY)
        cls.model_part.ProcessInfo[Kratos.DOMAIN_SIZE] = 3

        number_of_nodes = 10
        for id in range(1, number_of_nodes + 1):
            node = cls.model_part.CreateNewNode(id, id, id+1, id+2)
            node.SetSolutionStepValue(Kratos.VELOCITY, Kratos.Array3([id+3, id+4, id+5]))
            node.SetSolutionStepValue(Kratos.PRESSURE, id+3)
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

    def test_ContaienrDataAdd(self):
        a = KratosOA.HistoricalContainerVariableDataHolder(self.model_part)
        b = KratosOA.ElementPropertiesContainerVariableDataHolder(self.model_part)

        a.ReadDataFromContainerVariable(Kratos.VELOCITY)
        b.ReadDataFromContainerVariable(Kratos.PRESSURE)

        collective_1 = KratosOA.CollectiveVariableDataHolder()
        collective_1.AddVariableDataHolder(a)
        collective_1.AddVariableDataHolder(b)

        collective_2 = KratosOA.CollectiveVariableDataHolder()
        collective_2.AddVariableDataHolder(a)
        collective_2.AddVariableDataHolder(b)

        collective_3 = collective_1 * 10 + collective_2
        collective_3 += collective_1
        collective_4 = collective_3 + 4
        collective_4 += 10
        collective_4 += KratosOA.CollectiveVariableDataHolder(collective_1)

        collective_4.GetVariableDataHolders()[0].AssignDataToContainerVariable(Kratos.ACCELERATION)
        collective_4.GetVariableDataHolders()[1].AssignDataToContainerVariable(Kratos.DENSITY)

        for node in self.model_part.Nodes:
            self.assertVectorAlmostEqual(node.GetSolutionStepValue(Kratos.ACCELERATION), node.GetSolutionStepValue(Kratos.VELOCITY) * 13 + Kratos.Array3([14, 14, 14]), 12)
        for element in self.model_part.Elements:
            self.assertEqual(element.Properties[Kratos.DENSITY], element.Properties[Kratos.PRESSURE] * 13 + 14, 12)

    def test_ContaienrDataSub(self):
        a = KratosOA.HistoricalContainerVariableDataHolder(self.model_part)
        b = KratosOA.ElementPropertiesContainerVariableDataHolder(self.model_part)

        a.ReadDataFromContainerVariable(Kratos.VELOCITY)
        b.ReadDataFromContainerVariable(Kratos.PRESSURE)

        collective_1 = KratosOA.CollectiveVariableDataHolder()
        collective_1.AddVariableDataHolder(a)
        collective_1.AddVariableDataHolder(b)

        collective_2 = KratosOA.CollectiveVariableDataHolder()
        collective_2.AddVariableDataHolder(a)
        collective_2.AddVariableDataHolder(b)

        collective_3 = collective_1 * 10 - collective_2
        collective_3 -= collective_1
        collective_4 = collective_3 - 4
        collective_4 -= 10
        collective_4 -= KratosOA.CollectiveVariableDataHolder(collective_1)

        collective_4.GetVariableDataHolders()[0].AssignDataToContainerVariable(Kratos.ACCELERATION)
        collective_4.GetVariableDataHolders()[1].AssignDataToContainerVariable(Kratos.DENSITY)

        for node in self.model_part.Nodes:
            self.assertVectorAlmostEqual(node.GetSolutionStepValue(Kratos.ACCELERATION), node.GetSolutionStepValue(Kratos.VELOCITY) * 7 - Kratos.Array3([14, 14, 14]), 12)
        for element in self.model_part.Elements:
            self.assertEqual(element.Properties[Kratos.DENSITY], element.Properties[Kratos.PRESSURE] * 7 - 14, 12)

    def test_ContaienrDataMul(self):
        a = KratosOA.HistoricalContainerVariableDataHolder(self.model_part)
        b = KratosOA.ElementPropertiesContainerVariableDataHolder(self.model_part)

        a.ReadDataFromContainerVariable(Kratos.VELOCITY)
        b.ReadDataFromContainerVariable(Kratos.PRESSURE)

        collective_1 = KratosOA.CollectiveVariableDataHolder()
        collective_1.AddVariableDataHolder(a)
        collective_1.AddVariableDataHolder(b)

        collective_2 = KratosOA.CollectiveVariableDataHolder()
        collective_2.AddVariableDataHolder(a)
        collective_2.AddVariableDataHolder(b)

        collective_3 = collective_1 * 10
        collective_3 *= 2

        collective_3.GetVariableDataHolders()[0].AssignDataToContainerVariable(Kratos.ACCELERATION)
        collective_3.GetVariableDataHolders()[1].AssignDataToContainerVariable(Kratos.DENSITY)

        for node in self.model_part.Nodes:
            self.assertVectorAlmostEqual(node.GetSolutionStepValue(Kratos.ACCELERATION), node.GetSolutionStepValue(Kratos.VELOCITY) * 20, 12)
        for element in self.model_part.Elements:
            self.assertEqual(element.Properties[Kratos.DENSITY], element.Properties[Kratos.PRESSURE] * 20, 12)

    def test_ContaienrDataDiv(self):
        a = KratosOA.HistoricalContainerVariableDataHolder(self.model_part)
        b = KratosOA.ElementPropertiesContainerVariableDataHolder(self.model_part)

        a.ReadDataFromContainerVariable(Kratos.VELOCITY)
        b.ReadDataFromContainerVariable(Kratos.PRESSURE)

        collective_1 = KratosOA.CollectiveVariableDataHolder()
        collective_1.AddVariableDataHolder(a)
        collective_1.AddVariableDataHolder(b)

        collective_2 = KratosOA.CollectiveVariableDataHolder()
        collective_2.AddVariableDataHolder(a)
        collective_2.AddVariableDataHolder(b)

        collective_3 = collective_1 / 10
        collective_3 /= 2

        collective_3.GetVariableDataHolders()[0].AssignDataToContainerVariable(Kratos.ACCELERATION)
        collective_3.GetVariableDataHolders()[1].AssignDataToContainerVariable(Kratos.DENSITY)

        for node in self.model_part.Nodes:
            self.assertVectorAlmostEqual(node.GetSolutionStepValue(Kratos.ACCELERATION), node.GetSolutionStepValue(Kratos.VELOCITY) / 20, 12)
        for element in self.model_part.Elements:
            self.assertEqual(element.Properties[Kratos.DENSITY], element.Properties[Kratos.PRESSURE] / 20, 12)

    def test_ContaienrDataPow(self):
        a = KratosOA.HistoricalContainerVariableDataHolder(self.model_part)
        b = KratosOA.ElementPropertiesContainerVariableDataHolder(self.model_part)

        a.ReadDataFromContainerVariable(Kratos.VELOCITY)
        b.ReadDataFromContainerVariable(Kratos.PRESSURE)

        collective_1 = KratosOA.CollectiveVariableDataHolder()
        collective_1.AddVariableDataHolder(a)
        collective_1.AddVariableDataHolder(b)

        collective_2 = KratosOA.CollectiveVariableDataHolder()
        collective_2.AddVariableDataHolder(a)
        collective_2.AddVariableDataHolder(b)

        collective_3 = collective_1 ** 2
        collective_3 **= 2

        collective_3.GetVariableDataHolders()[0].AssignDataToContainerVariable(Kratos.ACCELERATION)
        collective_3.GetVariableDataHolders()[1].AssignDataToContainerVariable(Kratos.DENSITY)

        for node in self.model_part.Nodes:
            v = node.GetSolutionStepValue(Kratos.VELOCITY)
            self.assertVectorAlmostEqual(node.GetSolutionStepValue(Kratos.ACCELERATION), Kratos.Array3([v[0]**4, v[1]**4, v[2]**4]) , 12)
        for element in self.model_part.Elements:
            self.assertEqual(element.Properties[Kratos.DENSITY], element.Properties[Kratos.PRESSURE] ** 4, 12)

    def test_ContaienrDataNeg(self):
        a = KratosOA.HistoricalContainerVariableDataHolder(self.model_part)
        b = KratosOA.ElementPropertiesContainerVariableDataHolder(self.model_part)

        a.ReadDataFromContainerVariable(Kratos.VELOCITY)
        b.ReadDataFromContainerVariable(Kratos.PRESSURE)

        collective_1 = KratosOA.CollectiveVariableDataHolder()
        collective_1.AddVariableDataHolder(a)
        collective_1.AddVariableDataHolder(b)

        collective_2 = KratosOA.CollectiveVariableDataHolder()
        collective_2.AddVariableDataHolder(a)
        collective_2.AddVariableDataHolder(b)

        collective_3 = -collective_1

        collective_3.GetVariableDataHolders()[0].AssignDataToContainerVariable(Kratos.ACCELERATION)
        collective_3.GetVariableDataHolders()[1].AssignDataToContainerVariable(Kratos.DENSITY)

        for node in self.model_part.Nodes:
            v = node.GetSolutionStepValue(Kratos.VELOCITY)
            self.assertVectorAlmostEqual(node.GetSolutionStepValue(Kratos.ACCELERATION), v*(-1) , 12)
        for element in self.model_part.Elements:
            self.assertEqual(element.Properties[Kratos.DENSITY], -element.Properties[Kratos.PRESSURE], 12)

    def test_NormInf(self):
        a = KratosOA.HistoricalContainerVariableDataHolder(self.model_part)
        b = KratosOA.ElementPropertiesContainerVariableDataHolder(self.model_part)

        a.ReadDataFromContainerVariable(Kratos.VELOCITY)
        b.ReadDataFromContainerVariable(Kratos.PRESSURE)

        collective_1 = KratosOA.CollectiveVariableDataHolder([a, b])
        self.assertEqual(KratosOA.ContainerVariableDataHolderUtils.NormInf(collective_1), max(KratosOA.ContainerVariableDataHolderUtils.NormInf(a), KratosOA.ContainerVariableDataHolderUtils.NormInf(b)))

    def test_NormL2(self):
        a = KratosOA.HistoricalContainerVariableDataHolder(self.model_part)
        b = KratosOA.ElementPropertiesContainerVariableDataHolder(self.model_part)

        a.ReadDataFromContainerVariable(Kratos.VELOCITY)
        b.ReadDataFromContainerVariable(Kratos.PRESSURE)

        collective_1 = KratosOA.CollectiveVariableDataHolder([a, b])
        self.assertEqual(KratosOA.ContainerVariableDataHolderUtils.NormL2(collective_1), math.sqrt(KratosOA.ContainerVariableDataHolderUtils.NormL2(a)**2 + KratosOA.ContainerVariableDataHolderUtils.NormL2(b)**2))


    def test_InnerProduct(self):
        a = KratosOA.HistoricalContainerVariableDataHolder(self.model_part)
        b = KratosOA.ElementPropertiesContainerVariableDataHolder(self.model_part)

        collective_1 = KratosOA.CollectiveVariableDataHolder([a, b])

        a.ReadDataFromContainerVariable(Kratos.VELOCITY)
        b.ReadDataFromContainerVariable(Kratos.PRESSURE)

        self.assertEqual(KratosOA.ContainerVariableDataHolderUtils.InnerProduct(collective_1, collective_1), KratosOA.ContainerVariableDataHolderUtils.InnerProduct(a, a) + KratosOA.ContainerVariableDataHolderUtils.InnerProduct(b, b))

if __name__ == "__main__":
    Kratos.Tester.SetVerbosity(Kratos.Tester.Verbosity.PROGRESS)  # TESTS_OUTPUTS
    kratos_unittest.main()