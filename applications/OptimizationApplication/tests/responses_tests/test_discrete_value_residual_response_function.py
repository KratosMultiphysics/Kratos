from math import log
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
            "residual_type"          : "exact",
            "list_of_discrete_values": [0.0, 1.0, 2.0, 3.0]

        }""")
        cls.response_function = DiscreteValueResidualResponseFunction("geo_centroid", cls.model, parameters)
        cls.response_function.Initialize()
        cls.response_function.Check()
        cls.ref_value = cls.response_function.CalculateValue()

    def test_CalculateValue(self):
        value = 0.0
        for node in self.model_part.Nodes:
            v = node.GetValue(Kratos.PRESSURE)
            value += ((v) ** 2) * ((v - 1) ** 2) * ((v - 2) ** 2) * ((v - 3) ** 2)
        self.assertAlmostEqual(self.ref_value, value, 12)

    def test_CalculateGradient(self):
        ref_value = self.response_function.CalculateValue()
        analytical_gradient = Kratos.TensorAdaptors.VariableTensorAdaptor(self.model_part.Nodes, Kratos.PRESSURE)
        self.response_function.CalculateGradient({Kratos.PRESSURE: Kratos.TensorAdaptors.DoubleCombinedTensorAdaptor([analytical_gradient])})

        delta = 1e-9
        for i, node in enumerate(self.model_part.Nodes):
            node.SetValue(Kratos.PRESSURE, node.GetValue(Kratos.PRESSURE) + delta)
            fd_gradient = (self.response_function.CalculateValue() - ref_value) / delta
            node.SetValue(Kratos.PRESSURE, node.GetValue(Kratos.PRESSURE) - delta)
            self.assertAlmostEqual(fd_gradient, analytical_gradient.data[i], 5)

class TestDiscreteValueResidualResponseFunctionLogarithm(kratos_unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.model = Kratos.Model()
        cls.model_part = cls.model.CreateModelPart("test")
        with kratos_unittest.WorkFolderScope(".", __file__, True):
            ReadModelPart("../model_part_utils_test/quads", cls.model_part)

        for node in cls.model_part.Nodes:
            node.SetValue(Kratos.PRESSURE, (node.Id % 9 + 1) * 3 / 10)

        parameters = Kratos.Parameters("""{
            "evaluated_model_part_names"     : [
                "test"
            ],
            "container_type"         : "node_non_historical",
            "variable_name"          : "PRESSURE",
            "residual_type"          : "logarithm",
            "list_of_discrete_values": [0.0, 1.0, 2.0, 3.0]

        }""")
        cls.response_function = DiscreteValueResidualResponseFunction("geo_centroid", cls.model, parameters)
        cls.response_function.Initialize()
        cls.response_function.Check()
        cls.ref_value = cls.response_function.CalculateValue()

    def test_CalculateValue(self):
        value = 0.0
        for node in self.model_part.Nodes:
            v = node.GetValue(Kratos.PRESSURE)
            value += log((v) ** 2) + log((v - 1) ** 2) + log((v - 2) ** 2) + log((v - 3) ** 2)
        self.assertAlmostEqual(self.ref_value, value, 12)

    def test_CalculateGradient(self):
        ref_value = self.response_function.CalculateValue()
        analytical_gradient = Kratos.TensorAdaptors.VariableTensorAdaptor(self.model_part.Nodes, Kratos.PRESSURE)
        self.response_function.CalculateGradient({Kratos.PRESSURE: Kratos.TensorAdaptors.DoubleCombinedTensorAdaptor([analytical_gradient])})

        delta = 1e-8
        for i, node in enumerate(self.model_part.Nodes):
            node.SetValue(Kratos.PRESSURE, node.GetValue(Kratos.PRESSURE) + delta)
            fd_gradient = (self.response_function.CalculateValue() - ref_value) / delta
            node.SetValue(Kratos.PRESSURE, node.GetValue(Kratos.PRESSURE) - delta)
            self.assertAlmostEqual(fd_gradient, analytical_gradient.data[i], 5)

if __name__ == "__main__":
    kratos_unittest.main()