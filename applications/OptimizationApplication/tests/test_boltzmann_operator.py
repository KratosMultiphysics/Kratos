
import numpy
import KratosMultiphysics as Kratos
import KratosMultiphysics.OptimizationApplication as KratosOA

# Import KratosUnittest
import KratosMultiphysics.KratosUnittest as kratos_unittest

class TestBoltzmannOperator(kratos_unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.model = Kratos.Model()
        cls.model_part = cls.model.CreateModelPart("test")

        for i in range(10):
            node = cls.model_part.CreateNewNode(i + 1, 0, 0, 0)
            node.SetValue(Kratos.PRESSURE, i + 1 - 5)

        cls.nodal_ta = Kratos.TensorAdaptors.VariableTensorAdaptor(cls.model_part.Nodes, Kratos.PRESSURE)
        cls.nodal_ta.CollectData()

    def test_CalculateValue(self):
        max_boltzmann_operator = KratosOA.BoltzmannOperator(300)
        max_boltzmann_operator.Update(self.nodal_ta)
        self.assertAlmostEqual(numpy.max(self.nodal_ta.data), max_boltzmann_operator.CalculateValue())

        max_boltzmann_operator = KratosOA.BoltzmannOperator(-300)
        max_boltzmann_operator.Update(self.nodal_ta)
        self.assertAlmostEqual(numpy.min(self.nodal_ta.data), max_boltzmann_operator.CalculateValue())

    def test_CalculateGradient(self):
        boltzmann_operator = KratosOA.BoltzmannOperator(3.0)
        boltzmann_operator.Update(self.nodal_ta)
        ref_value = boltzmann_operator.CalculateValue()
        analytical_gradient = boltzmann_operator.CalculateGradient()

        delta = 1e-9
        for i in range(self.nodal_ta.data.shape[0]):
            self.nodal_ta.data[i] += delta
            boltzmann_operator.Update(self.nodal_ta)
            sensitivity = (boltzmann_operator.CalculateValue() - ref_value) / delta
            self.nodal_ta.data[i] -= delta
            self.assertAlmostEqual(sensitivity, analytical_gradient.data[i], 5)

if __name__ == "__main__":
    kratos_unittest.main()