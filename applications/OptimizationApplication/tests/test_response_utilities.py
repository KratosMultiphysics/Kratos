
from math import log

import KratosMultiphysics as Kratos
import KratosMultiphysics.OptimizationApplication as KratosOA
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem import OptimizationProblem
from KratosMultiphysics.OptimizationApplication.utilities.component_data_view import ComponentDataView
from KratosMultiphysics.OptimizationApplication.responses.discrete_value_residual_response_function import DiscreteValueResidualResponseFunction
from KratosMultiphysics.OptimizationApplication.utilities.response_utilities import EvaluateResponseExpression

# Import KratosUnittest
import KratosMultiphysics.KratosUnittest as kratos_unittest

class TestResponseUtilities(kratos_unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.model = Kratos.Model()
        cls.model_part = cls.model.CreateModelPart("test")

        for i in range(10):
            node: Kratos.Node = cls.model_part.CreateNewNode(i + 1, i, i + 1, i + 2)
            node.SetValue(Kratos.PRESSURE, (i + 1) / 100.0)

        r1_params = Kratos.Parameters("""{
            "evaluated_model_part_names"     : [
                "test"
            ],
            "container_type"         : "node_non_historical",
            "variable_name"          : "PRESSURE",
            "residual_type"          : "exact",
            "list_of_discrete_values": [-1.0e-2, -2.0e-2, -3.0e-2]
        }""")
        cls.r1 = DiscreteValueResidualResponseFunction("r1", cls.model, r1_params)
        cls.r1.Initialize()

        r2_params = Kratos.Parameters("""{
            "evaluated_model_part_names"     : [
                "test"
            ],
            "container_type"         : "node_non_historical",
            "variable_name"          : "PRESSURE",
            "residual_type"          : "exact",
            "list_of_discrete_values": [-1e-3, -2e-3, -3e-3]
        }""")
        cls.r2 = DiscreteValueResidualResponseFunction("r2", cls.model, r2_params)
        cls.r2.Initialize()

    def setUp(self) -> None:
        self.optimization_problem = OptimizationProblem()
        self.optimization_problem.AddComponent(self.r1)
        ComponentDataView(self.r1, self.optimization_problem).SetDataBuffer(1)
        self.optimization_problem.AddComponent(self.r2)
        ComponentDataView(self.r2, self.optimization_problem).SetDataBuffer(1)

    def test_LiteralResponseCalculateValue1(self):
        eval_resp = EvaluateResponseExpression(self.model, "4.0 + 6.0", self.optimization_problem)
        eval_resp.Initialize()
        eval_resp.Initialize()
        eval_value = eval_resp.CalculateValue()
        self.assertEqual(eval_value, 10.0)

        self.assertFalse(self.optimization_problem.HasResponse(eval_resp))
        self.assertFalse(self.optimization_problem.HasResponse("4.0"))
        self.assertFalse(self.optimization_problem.HasResponse("6.0"))
        self.assertFalse(self.optimization_problem.HasResponse("(4.0+6.0)"))

    def test_LiteralResponseCalculateValue2(self):
        eval_resp = EvaluateResponseExpression(self.model, "4.0 * 6.0", self.optimization_problem)
        eval_resp.Initialize()
        eval_value = eval_resp.CalculateValue()
        self.assertEqual(eval_value, 24.0)

        self.assertFalse(self.optimization_problem.HasResponse(eval_resp))
        self.assertFalse(self.optimization_problem.HasResponse("4.0"))
        self.assertFalse(self.optimization_problem.HasResponse("6.0"))
        self.assertFalse(self.optimization_problem.HasResponse("(4.0*6.0)"))

    def test_LiteralResponseCalculateValue3(self):
        eval_resp = EvaluateResponseExpression(self.model, "4.0 / 8.0", self.optimization_problem)
        eval_resp.Initialize()
        eval_value = eval_resp.CalculateValue()
        self.assertEqual(eval_value, 0.5)

        self.assertFalse(self.optimization_problem.HasResponse(eval_resp))
        self.assertFalse(self.optimization_problem.HasResponse("4.0"))
        self.assertFalse(self.optimization_problem.HasResponse("8.0"))
        self.assertFalse(self.optimization_problem.HasResponse("(4.0/6.0)"))

    def test_LiteralResponseCalculateValue4(self):
        eval_resp = EvaluateResponseExpression(self.model, "4.0 - 8.0", self.optimization_problem)
        eval_resp.Initialize()
        eval_value = eval_resp.CalculateValue()
        self.assertEqual(eval_value, -4.0)

        self.assertFalse(self.optimization_problem.HasResponse(eval_resp))
        self.assertFalse(self.optimization_problem.HasResponse("4.0"))
        self.assertFalse(self.optimization_problem.HasResponse("8.0"))
        self.assertFalse(self.optimization_problem.HasResponse("(4.0-6.0)"))

    def test_LiteralResponseCalculateValue5(self):
        eval_resp = EvaluateResponseExpression(self.model, "4.0 ^ 2.0", self.optimization_problem)
        eval_resp.Initialize()
        eval_value = eval_resp.CalculateValue()
        self.assertEqual(eval_value, 16.0)

        self.assertFalse(self.optimization_problem.HasResponse(eval_resp))
        self.assertFalse(self.optimization_problem.HasResponse("4.0"))
        self.assertFalse(self.optimization_problem.HasResponse("8.0"))
        self.assertFalse(self.optimization_problem.HasResponse("(4.0^6.0)"))

    def test_LiteralResponseCalculateValue6(self):
        eval_resp = EvaluateResponseExpression(self.model, "4.0 - 8.0 + 3.0 * 2.0 / 5.0 ^ 2", self.optimization_problem)
        eval_resp.Initialize()
        eval_value = eval_resp.CalculateValue()
        self.assertEqual(eval_value, -3.76)

        self.assertFalse(self.optimization_problem.HasResponse(eval_resp))
        self.assertFalse(self.optimization_problem.HasResponse("4.0"))
        self.assertFalse(self.optimization_problem.HasResponse("8.0"))
        self.assertFalse(self.optimization_problem.HasResponse("3.0"))
        self.assertFalse(self.optimization_problem.HasResponse("2.0"))
        self.assertFalse(self.optimization_problem.HasResponse("5.0"))

    def test_CombinedResponseCalculateValue1(self):
        eval_resp = EvaluateResponseExpression(self.model, "r1", self.optimization_problem)
        eval_resp.Initialize()
        eval_value = eval_resp.CalculateValue()
        self.assertEqual(eval_value, self.r1.CalculateValue())

        self.assertTrue(self.optimization_problem.HasResponse(eval_resp))
        self.assertTrue(self.optimization_problem.HasResponse("r1"))
        self.assertEqual(eval_resp.GetName(), "r1")

        # following is the resultant response function, hence the value storage is managed by the ResponseRoutine
        self.assertFalse(ComponentDataView(self.optimization_problem.GetResponse("r1"), self.optimization_problem).GetBufferedData().HasValue("value"))

    def test_CombinedResponseCalculateValue2(self):
        eval_resp = EvaluateResponseExpression(self.model, "r1 + r2 + r2", self.optimization_problem)
        eval_resp.Initialize()
        eval_value = eval_resp.CalculateValue()
        self.assertAlmostEqual(eval_value, self.r1.CalculateValue() + 2.0 * self.r2.CalculateValue(), 12)

        # followings are intermediate responses, or leaf responses. Therefore they need to be
        # in the optimization problem
        self.assertTrue(self.optimization_problem.HasResponse("r1"))
        self.assertTrue(self.optimization_problem.HasResponse("r2"))
        self.assertTrue(self.optimization_problem.HasResponse("(r1+r2)"))

        # followings are intermediate responses, so Algorithm, ResponseRoutine do not see them. Therefore
        # the value storage is managed by the ResponseFunction it self.
        self.assertTrue(ComponentDataView(self.optimization_problem.GetResponse("r1"), self.optimization_problem).GetBufferedData().HasValue("value"))
        self.assertTrue(ComponentDataView(self.optimization_problem.GetResponse("r2"), self.optimization_problem).GetBufferedData().HasValue("value"))
        self.assertTrue(ComponentDataView(self.optimization_problem.GetResponse("(r1+r2)"), self.optimization_problem).GetBufferedData().HasValue("value"))

        # following is the resultant response function, hence the value storage is managed by the ResponseRoutine
        self.assertFalse(self.optimization_problem.HasResponse(eval_resp))
        self.assertFalse(self.optimization_problem.HasResponse("((r1+r2)+r2)"))

    def test_CombinedResponseCalculateValue3(self):
        eval_resp = EvaluateResponseExpression(self.model, "r1 + r2 + 4.0", self.optimization_problem)
        eval_resp.Initialize()
        eval_value = eval_resp.CalculateValue()
        self.assertEqual(eval_value, self.r1.CalculateValue() + self.r2.CalculateValue() + 4.0)

        # followings are intermediate responses, or leaf responses. Therefore they need to be
        # in the optimization problem
        self.assertTrue(self.optimization_problem.HasResponse("r1"))
        self.assertTrue(self.optimization_problem.HasResponse("r2"))
        self.assertTrue(self.optimization_problem.HasResponse("(r1+r2)"))
        self.assertFalse(self.optimization_problem.HasResponse("4.0"))

        # followings are intermediate responses, so Algorithm, ResponseRoutine do not see them. Therefore
        # the value storage is managed by the ResponseFunction it self.
        self.assertTrue(ComponentDataView(self.optimization_problem.GetResponse("r1"), self.optimization_problem).GetBufferedData().HasValue("value"))
        self.assertTrue(ComponentDataView(self.optimization_problem.GetResponse("r2"), self.optimization_problem).GetBufferedData().HasValue("value"))
        self.assertTrue(ComponentDataView(self.optimization_problem.GetResponse("(r1+r2)"), self.optimization_problem).GetBufferedData().HasValue("value"))

        # following is the resultant response function, hence the value storage is managed by the ResponseRoutine
        self.assertFalse(self.optimization_problem.HasResponse(eval_resp))
        self.assertFalse(self.optimization_problem.HasResponse("((r1+r2)+4.0)"))

    def test_CombinedResponseCalculateValue4(self):
        eval_resp = EvaluateResponseExpression(self.model, "3.5 + r1 * 2.0 + r2 / 3.0 - 4.0", self.optimization_problem)
        eval_resp.Initialize()
        eval_value = eval_resp.CalculateValue()
        self.assertEqual(eval_value, self.r1.CalculateValue() * 2.0 + self.r2.CalculateValue() / 3.0 - 4.0 + 3.5)

        # followings are intermediate responses, or leaf responses. Therefore they need to be
        # in the optimization problem
        self.assertTrue(self.optimization_problem.HasResponse("r1"))
        self.assertTrue(self.optimization_problem.HasResponse("r2"))
        self.assertTrue(self.optimization_problem.HasResponse("(r1*2.0)"))
        self.assertTrue(self.optimization_problem.HasResponse("(r2÷3.0)"))
        self.assertTrue(self.optimization_problem.HasResponse("(3.5+(r1*2.0))"))
        self.assertTrue(self.optimization_problem.HasResponse("((3.5+(r1*2.0))+(r2÷3.0))"))

        self.assertTrue(ComponentDataView(self.optimization_problem.GetResponse("r1"), self.optimization_problem).GetBufferedData().HasValue("value"))
        self.assertTrue(ComponentDataView(self.optimization_problem.GetResponse("r2"), self.optimization_problem).GetBufferedData().HasValue("value"))
        self.assertTrue(ComponentDataView(self.optimization_problem.GetResponse("(r1*2.0)"), self.optimization_problem).GetBufferedData().HasValue("value"))
        self.assertTrue(ComponentDataView(self.optimization_problem.GetResponse("(r2÷3.0)"), self.optimization_problem).GetBufferedData().HasValue("value"))
        self.assertTrue(ComponentDataView(self.optimization_problem.GetResponse("(3.5+(r1*2.0))"), self.optimization_problem).GetBufferedData().HasValue("value"))
        self.assertTrue(ComponentDataView(self.optimization_problem.GetResponse("((3.5+(r1*2.0))+(r2÷3.0))"), self.optimization_problem).GetBufferedData().HasValue("value"))

        # following is the resultant response function, hence the value storage is managed by the ResponseRoutine
        self.assertFalse(self.optimization_problem.HasResponse(eval_resp))
        self.assertFalse(self.optimization_problem.HasResponse("(((3.5+(r1*2.0))+(r2/3.0))-4.0)"))

    def test_CombinedResponseCalculateValue5(self):
        eval_resp = EvaluateResponseExpression(self.model, "3.5 + r1 * 2.0 * r2 / 3.0 - 4.0 + r1 ^ r2", self.optimization_problem)
        eval_resp.Initialize()
        eval_value = eval_resp.CalculateValue()
        self.assertEqual(eval_value, self.r1.CalculateValue() * 2.0 * self.r2.CalculateValue() / 3.0 - 4.0 + 3.5 + self.r1.CalculateValue() ** self.r2.CalculateValue())

        # followings are intermediate responses, or leaf responses. Therefore they need to be
        # in the optimization problem
        self.assertTrue(self.optimization_problem.HasResponse("r1"))
        self.assertTrue(self.optimization_problem.HasResponse("r2"))
        self.assertTrue(self.optimization_problem.HasResponse("(r1*2.0)"))
        self.assertTrue(self.optimization_problem.HasResponse("((r1*2.0)*r2)"))
        self.assertTrue(self.optimization_problem.HasResponse("(((r1*2.0)*r2)÷3.0)"))
        self.assertTrue(self.optimization_problem.HasResponse("(r1^r2)"))
        self.assertTrue(self.optimization_problem.HasResponse("(3.5+(((r1*2.0)*r2)÷3.0))"))
        self.assertTrue(self.optimization_problem.HasResponse("((3.5+(((r1*2.0)*r2)÷3.0))-4.0)"))

        self.assertTrue(ComponentDataView(self.optimization_problem.GetResponse("r1"), self.optimization_problem).GetBufferedData().HasValue("value"))
        self.assertTrue(ComponentDataView(self.optimization_problem.GetResponse("r2"), self.optimization_problem).GetBufferedData().HasValue("value"))
        self.assertTrue(ComponentDataView(self.optimization_problem.GetResponse("(r1*2.0)"), self.optimization_problem).GetBufferedData().HasValue("value"))
        self.assertTrue(ComponentDataView(self.optimization_problem.GetResponse("((r1*2.0)*r2)"), self.optimization_problem).GetBufferedData().HasValue("value"))
        self.assertTrue(ComponentDataView(self.optimization_problem.GetResponse("(((r1*2.0)*r2)÷3.0)"), self.optimization_problem).GetBufferedData().HasValue("value"))
        self.assertTrue(ComponentDataView(self.optimization_problem.GetResponse("(r1^r2)"), self.optimization_problem).GetBufferedData().HasValue("value"))
        self.assertTrue(ComponentDataView(self.optimization_problem.GetResponse("(3.5+(((r1*2.0)*r2)÷3.0))"), self.optimization_problem).GetBufferedData().HasValue("value"))
        self.assertTrue(ComponentDataView(self.optimization_problem.GetResponse("((3.5+(((r1*2.0)*r2)÷3.0))-4.0)"), self.optimization_problem).GetBufferedData().HasValue("value"))

        # following is the resultant response function, hence the value storage is managed by the ResponseRoutine
        self.assertFalse(self.optimization_problem.HasResponse(eval_resp))
        self.assertFalse(self.optimization_problem.HasResponse(eval_resp.GetName()))

    def test_BracketResponseCalculateValue1(self):
        eval_resp = EvaluateResponseExpression(self.model, "(4.0)", self.optimization_problem)
        eval_resp.Initialize()
        eval_resp.Initialize()
        eval_value = eval_resp.CalculateValue()
        self.assertEqual(eval_value, 4.0)

    def test_BracketResponseCalculateValue2(self):
        eval_resp = EvaluateResponseExpression(self.model, "(4.0 + (6.0 / 3.0 + 3 * (2 + 8))) * 2.0", self.optimization_problem)
        eval_resp.Initialize()
        eval_value = eval_resp.CalculateValue()
        self.assertEqual(eval_value, 72.0)

    def test_BracketResponseCalculateValue3(self):
        eval_resp = EvaluateResponseExpression(self.model, "(4.0 + (r1 * 2) * (r2 ^ 2) + log((2+6)*(3+(4-2)/4)+6*r1) + (6.0 / 3.0 + 3 * (2 + 8))) * 2.0", self.optimization_problem)
        eval_resp.Initialize()
        eval_value = eval_resp.CalculateValue()
        self.assertEqual(eval_value, (4.0 + (self.r1.CalculateValue() * 2) * (self.r2.CalculateValue() ** 2) + log((2+6)*(3+(4-2)/4)+6*self.r1.CalculateValue()) + (6.0 / 3.0 + 3 * (2 + 8))) * 2.0)

        # followings are intermediate responses, or leaf responses. Therefore they need to be
        # in the optimization problem
        self.assertTrue(self.optimization_problem.HasResponse("r1"))
        self.assertTrue(self.optimization_problem.HasResponse("r2"))
        self.assertTrue(self.optimization_problem.HasResponse("(r1*2.0)"))
        self.assertTrue(self.optimization_problem.HasResponse("(r2^2.0)"))
        self.assertTrue(self.optimization_problem.HasResponse("((r1*2.0)*(r2^2.0))"))
        self.assertTrue(self.optimization_problem.HasResponse("(6.0*r1)"))
        self.assertTrue(self.optimization_problem.HasResponse("((4.0+((r1*2.0)*(r2^2.0)))+log((((2.0+6.0)*(3.0+((4.0-2.0)÷4.0)))+(6.0*r1))))"))
        self.assertTrue(self.optimization_problem.HasResponse("(((4.0+((r1*2.0)*(r2^2.0)))+log((((2.0+6.0)*(3.0+((4.0-2.0)÷4.0)))+(6.0*r1))))+((6.0÷3.0)+(3.0*(2.0+8.0))))"))
        self.assertTrue(self.optimization_problem.HasResponse("log((((2.0+6.0)*(3.0+((4.0-2.0)÷4.0)))+(6.0*r1)))"))

        self.assertTrue(ComponentDataView(self.optimization_problem.GetResponse("r1"), self.optimization_problem).GetBufferedData().HasValue("value"))
        self.assertTrue(ComponentDataView(self.optimization_problem.GetResponse("r2"), self.optimization_problem).GetBufferedData().HasValue("value"))
        self.assertTrue(ComponentDataView(self.optimization_problem.GetResponse("(r1*2.0)"), self.optimization_problem).GetBufferedData().HasValue("value"))
        self.assertTrue(ComponentDataView(self.optimization_problem.GetResponse("(r2^2.0)"), self.optimization_problem).GetBufferedData().HasValue("value"))
        self.assertTrue(ComponentDataView(self.optimization_problem.GetResponse("((r1*2.0)*(r2^2.0))"), self.optimization_problem).GetBufferedData().HasValue("value"))
        self.assertTrue(ComponentDataView(self.optimization_problem.GetResponse("(6.0*r1)"), self.optimization_problem).GetBufferedData().HasValue("value"))
        self.assertTrue(ComponentDataView(self.optimization_problem.GetResponse("((4.0+((r1*2.0)*(r2^2.0)))+log((((2.0+6.0)*(3.0+((4.0-2.0)÷4.0)))+(6.0*r1))))"), self.optimization_problem).GetBufferedData().HasValue("value"))
        self.assertTrue(ComponentDataView(self.optimization_problem.GetResponse("(((4.0+((r1*2.0)*(r2^2.0)))+log((((2.0+6.0)*(3.0+((4.0-2.0)÷4.0)))+(6.0*r1))))+((6.0÷3.0)+(3.0*(2.0+8.0))))"), self.optimization_problem).GetBufferedData().HasValue("value"))
        self.assertTrue(ComponentDataView(self.optimization_problem.GetResponse("log((((2.0+6.0)*(3.0+((4.0-2.0)÷4.0)))+(6.0*r1)))"), self.optimization_problem).GetBufferedData().HasValue("value"))

        # following is the resultant response function, hence the value storage is managed by the ResponseRoutine
        self.assertFalse(self.optimization_problem.HasResponse(eval_resp))
        self.assertFalse(self.optimization_problem.HasResponse(eval_resp.GetName()))

if __name__ == "__main__":
    kratos_unittest.main()