import math

import KratosMultiphysics as Kratos
import KratosMultiphysics.OptimizationApplication as KratosOA

# Import KratosUnittest
import KratosMultiphysics.KratosUnittest as kratos_unittest


class TestContainerDataUtils(kratos_unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.model = Kratos.Model()
        cls.model_part = cls.model.CreateModelPart("test")
        cls.model_part.AddNodalSolutionStepVariable(Kratos.PRESSURE)
        cls.model_part.AddNodalSolutionStepVariable(Kratos.DENSITY)
        cls.model_part.AddNodalSolutionStepVariable(Kratos.VELOCITY)
        cls.model_part.ProcessInfo[Kratos.DOMAIN_SIZE] = 3

        number_of_nodes = 10
        for id in range(1, number_of_nodes + 1):
            node = cls.model_part.CreateNewNode(id, id, id+1, id+2)
            node.SetSolutionStepValue(Kratos.VELOCITY, Kratos.Array3([id+3, id+4, id+5]))
            node.SetSolutionStepValue(Kratos.PRESSURE, id+3)
            node.SetSolutionStepValue(Kratos.DENSITY, id+4)

    def test_NormInf(self):
        a = KratosOA.HistoricalContainerData(self.model_part)

        a.ReadDataFromContainerVariable(Kratos.PRESSURE)
        self.assertEqual(KratosOA.ContainerDataUtils.NormInf(a), 13)

        a.ReadDataFromContainerVariable(Kratos.VELOCITY)
        self.assertEqual(KratosOA.ContainerDataUtils.NormInf(a), 15)

    def test_EntityMaxNormL2(self):
        a = KratosOA.HistoricalContainerData(self.model_part)

        a.ReadDataFromContainerVariable(Kratos.PRESSURE)
        self.assertEqual(KratosOA.ContainerDataUtils.EntityMaxNormL2(a), 13)

        a.ReadDataFromContainerVariable(Kratos.VELOCITY)
        self.assertEqual(KratosOA.ContainerDataUtils.EntityMaxNormL2(a), math.sqrt(15**2 + 14**2 + 13**2))

    def test_InnerProduct(self):
        a = KratosOA.HistoricalContainerData(self.model_part)
        b = KratosOA.HistoricalContainerData(self.model_part)

        a.ReadDataFromContainerVariable(Kratos.PRESSURE)
        b.ReadDataFromContainerVariable(Kratos.DENSITY)

        self.assertEqual(KratosOA.ContainerDataUtils.InnerProduct(a, b), 890)

if __name__ == "__main__":
    Kratos.Tester.SetVerbosity(Kratos.Tester.Verbosity.PROGRESS)  # TESTS_OUTPUTS
    kratos_unittest.main()