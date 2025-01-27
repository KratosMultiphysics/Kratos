
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

        cls.boltzmann_operator = KratosOA.NodalBoltzmannOperator(3.0, 12.0)

        nodal_exp = Kratos.Expression.NodalExpression(cls.model_part)
        Kratos.Expression.VariableExpressionIO.Read(nodal_exp, Kratos.PRESSURE, False)
        cls.numpy_exp = nodal_exp.Evaluate()

        cls.c_exp = nodal_exp.Clone()
        Kratos.Expression.CArrayExpressionIO.Move(cls.c_exp, cls.numpy_exp)
        cls.boltzmann_operator.Update(cls.c_exp)

    def test_CalculateValue(self):
        max_boltzmann_operator = KratosOA.NodalBoltzmannOperator(300, 12)
        max_boltzmann_operator.Update(self.c_exp)
        self.assertAlmostEqual(5.0, max_boltzmann_operator.CalculateValue())

        max_boltzmann_operator = KratosOA.NodalBoltzmannOperator(-300, 12)
        max_boltzmann_operator.Update(self.c_exp)
        self.assertAlmostEqual(-4.0, max_boltzmann_operator.CalculateValue())

    def test_CalculateGradient(self):
        analytical_gradient = self.boltzmann_operator.CalculateGradient().Evaluate()
        ref_value = self.boltzmann_operator.CalculateValue()

        delta = 1e-9
        for i in range(self.numpy_exp.shape[0]):
            self.numpy_exp[i] += delta
            self.boltzmann_operator.Update(self.c_exp)
            sensitivity = (self.boltzmann_operator.CalculateValue() - ref_value) / delta
            self.numpy_exp[i] -= delta
            self.assertAlmostEqual(sensitivity, analytical_gradient[i], 5)

if __name__ == "__main__":
    kratos_unittest.main()