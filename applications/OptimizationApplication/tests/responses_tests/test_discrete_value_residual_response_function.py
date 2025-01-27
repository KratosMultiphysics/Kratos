from math import exp
import KratosMultiphysics as Kratos
import KratosMultiphysics.OptimizationApplication as KratosOA

import KratosMultiphysics.KratosUnittest as kratos_unittest
from KratosMultiphysics.testing.utilities import ReadModelPart
from KratosMultiphysics.OptimizationApplication.responses.discrete_value_residual_response_function import DiscreteValueResidualResponseFunction

class TestDiscreteValueResidualResponseFunctionExact(kratos_unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.model = Kratos.Model()
        cls.model_part = cls.model.CreateModelPart("test")
        with kratos_unittest.WorkFolderScope(".", __file__, True):
            ReadModelPart("../model_part_utils_test/quads", cls.model_part)

        for node in cls.model_part.Nodes:
            node.SetValue(Kratos.PRESSURE, (node.Id % 11) * 3 / 10)

        parameters = Kratos.Parameters("""{
            "evaluated_model_part_names"     : [
                "test"
            ],
            "container_type"         : "node_non_historical",
            "variable_name"          : "PRESSURE",
            "beta_coefficient"       : 2.0,
            "scaling_factor"         : 3.0,
            "list_of_discrete_values": [0.0, 1.0, 2.0, 3.0]

        }""")
        cls.response_function = DiscreteValueResidualResponseFunction("geo_centroid", cls.model, parameters)
        cls.response_function.Initialize()
        cls.response_function.Check()
        cls.ref_value = cls.response_function.CalculateValue()

    def test_CalculateValue(self):
        numerator = 0.0
        denominator = 0.0
        for node in self.model_part.Nodes:
            v = node.GetValue(Kratos.PRESSURE)
            v_res = ((v) ** 2) * ((v - 1) ** 2) * ((v - 2) ** 2) * ((v - 3) ** 2) / 3.0
            numerator += v_res * exp(2.0 * v_res)
            denominator += exp(2.0 * v_res)
        self.assertAlmostEqual(self.ref_value, 3.0 * numerator / denominator, 10)

    def test_CalculateGradient(self):
        ref_value = self.response_function.CalculateValue()
        analytical_gradient = Kratos.Expression.NodalExpression(self.model_part)
        self.response_function.CalculateGradient({Kratos.PRESSURE: KratosOA.CollectiveExpression([analytical_gradient])})
        analytical_gradient = analytical_gradient.Evaluate()

        delta = 1e-9
        for i, node in enumerate(self.model_part.Nodes):
            node.SetValue(Kratos.PRESSURE, node.GetValue(Kratos.PRESSURE) + delta)
            fd_gradient = (self.response_function.CalculateValue() - ref_value) / delta
            node.SetValue(Kratos.PRESSURE, node.GetValue(Kratos.PRESSURE) - delta)
            self.assertAlmostEqual(fd_gradient, analytical_gradient[i], 5)

if __name__ == "__main__":
    kratos_unittest.main()