import numpy as np
from math import log

import KratosMultiphysics as Kratos
import KratosMultiphysics.OptimizationApplication as KratosOA
from KratosMultiphysics.OptimizationApplication.responses.response_function import ResponseFunction
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem import OptimizationProblem
from KratosMultiphysics.OptimizationApplication.utilities.component_data_view import ComponentDataView
from KratosMultiphysics.OptimizationApplication.responses.discrete_value_residual_response_function import DiscreteValueResidualResponseFunction
from KratosMultiphysics.OptimizationApplication.utilities.response_utilities import EvaluateResponseExpression
from KratosMultiphysics.OptimizationApplication.utilities.union_utilities import SupportedSensitivityFieldVariableTypes

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
            node.SetValue(Kratos.TEMPERATURE, (i + 10) / 100.0)

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

        r3_params = Kratos.Parameters("""{
            "evaluated_model_part_names"     : [
                "test"
            ],
            "container_type"         : "node_non_historical",
            "variable_name"          : "TEMPERATURE",
            "residual_type"          : "exact",
            "list_of_discrete_values": [-2e-3, -3e-3, -4e-3]
        }""")
        cls.r3 = DiscreteValueResidualResponseFunction("r3", cls.model, r3_params)
        cls.r3.Initialize()

    def __CheckGradients(self, sensitivity_variables: 'list[SupportedSensitivityFieldVariableTypes]', response_function: ResponseFunction, analytical_lambda, delta: float, precision: int) -> None:
        ref_value = analytical_lambda()
        self.assertAlmostEqual(ref_value, response_function.CalculateValue(), precision)

        physical_variable_collective_expressions: 'dict[SupportedSensitivityFieldVariableTypes, KratosOA.CollectiveExpression]' = {}
        fd_gradients: 'dict[SupportedSensitivityFieldVariableTypes, list[float]]' = {}
        for sensitivity_variable in sensitivity_variables:
            nodal_exp = Kratos.Expression.NodalExpression(self.model_part)
            Kratos.Expression.LiteralExpressionIO.SetData(nodal_exp, 0.0)
            physical_variable_collective_expressions[sensitivity_variable] = KratosOA.CollectiveExpression([nodal_exp])
            fd_gradients[sensitivity_variable] = []

        response_function.CalculateGradient(physical_variable_collective_expressions)

        for node in self.model_part.Nodes:
            for sensitivity_variable in sensitivity_variables:
                node.SetValue(sensitivity_variable, node.GetValue(sensitivity_variable) + delta)
                fd_gradients[sensitivity_variable].append((analytical_lambda() - ref_value) / delta)
                node.SetValue(sensitivity_variable, node.GetValue(sensitivity_variable) - delta)

        for sensitivity_variable in sensitivity_variables:
            nodal_exp = Kratos.Expression.NodalExpression(self.model_part)
            Kratos.Expression.CArrayExpressionIO.Read(nodal_exp, np.array(fd_gradients[sensitivity_variable], np.float64))
            self.assertAlmostEqual(Kratos.Expression.Utils.NormInf(nodal_exp - physical_variable_collective_expressions[sensitivity_variable].GetContainerExpressions()[0]), 0.0, precision)

    def setUp(self) -> None:
        self.optimization_problem = OptimizationProblem()
        self.optimization_problem.AddComponent(self.r1)
        ComponentDataView(self.r1, self.optimization_problem).SetDataBuffer(1)
        self.optimization_problem.AddComponent(self.r2)
        ComponentDataView(self.r2, self.optimization_problem).SetDataBuffer(1)
        self.optimization_problem.AddComponent(self.r3)
        ComponentDataView(self.r3, self.optimization_problem).SetDataBuffer(1)

    def test_LiteralResponseCalculateValue1(self):
        eval_resp = EvaluateResponseExpression("4.0 + 6.0", self.optimization_problem)
        eval_resp.Initialize()
        eval_value = eval_resp.CalculateValue()
        self.assertEqual(eval_value, 10.0)

        self.assertFalse(self.optimization_problem.HasResponse(eval_resp))
        self.assertFalse(self.optimization_problem.HasResponse("4.0"))
        self.assertFalse(self.optimization_problem.HasResponse("6.0"))
        self.assertFalse(self.optimization_problem.HasResponse("(4.0+6.0)"))

        self.__CheckGradients([Kratos.PRESSURE, Kratos.TEMPERATURE], eval_resp, lambda : 4.0 + 6.0, 1e-8, 12)

    def test_LiteralResponseCalculateValue2(self):
        eval_resp = EvaluateResponseExpression("4.0 * 6.0", self.optimization_problem)
        eval_resp.Initialize()
        eval_value = eval_resp.CalculateValue()
        self.assertEqual(eval_value, 24.0)

        self.assertFalse(self.optimization_problem.HasResponse(eval_resp))
        self.assertFalse(self.optimization_problem.HasResponse("4.0"))
        self.assertFalse(self.optimization_problem.HasResponse("6.0"))
        self.assertFalse(self.optimization_problem.HasResponse("(4.0*6.0)"))

        self.__CheckGradients([Kratos.PRESSURE, Kratos.TEMPERATURE], eval_resp, lambda : 4.0 * 6.0, 1e-8, 12)

    def test_LiteralResponseCalculateValue3(self):
        eval_resp = EvaluateResponseExpression("4.0 / 8.0", self.optimization_problem)
        eval_resp.Initialize()
        eval_value = eval_resp.CalculateValue()
        self.assertEqual(eval_value, 0.5)

        self.assertFalse(self.optimization_problem.HasResponse(eval_resp))
        self.assertFalse(self.optimization_problem.HasResponse("4.0"))
        self.assertFalse(self.optimization_problem.HasResponse("8.0"))
        self.assertFalse(self.optimization_problem.HasResponse("(4.0/8.0)"))

        self.__CheckGradients([Kratos.PRESSURE, Kratos.TEMPERATURE], eval_resp, lambda : 4.0 / 8.0, 1e-8, 12)

    def test_LiteralResponseCalculateValue4(self):
        eval_resp = EvaluateResponseExpression("4.0 - 8.0", self.optimization_problem)
        eval_resp.Initialize()
        eval_value = eval_resp.CalculateValue()
        self.assertEqual(eval_value, -4.0)

        self.assertFalse(self.optimization_problem.HasResponse(eval_resp))
        self.assertFalse(self.optimization_problem.HasResponse("4.0"))
        self.assertFalse(self.optimization_problem.HasResponse("8.0"))
        self.assertFalse(self.optimization_problem.HasResponse("(4.0-8.0)"))

        self.__CheckGradients([Kratos.PRESSURE, Kratos.TEMPERATURE], eval_resp, lambda : 4.0 - 8.0, 1e-8, 12)

    def test_LiteralResponseCalculateValue5(self):
        eval_resp = EvaluateResponseExpression("4.0 ^ 2.0", self.optimization_problem)
        eval_resp.Initialize()
        eval_value = eval_resp.CalculateValue()
        self.assertEqual(eval_value, 16.0)

        self.assertFalse(self.optimization_problem.HasResponse(eval_resp))
        self.assertFalse(self.optimization_problem.HasResponse("4.0"))
        self.assertFalse(self.optimization_problem.HasResponse("8.0"))
        self.assertFalse(self.optimization_problem.HasResponse("(4.0^2.0)"))

        self.__CheckGradients([Kratos.PRESSURE, Kratos.TEMPERATURE], eval_resp, lambda : 4.0 ** 2.0, 1e-8, 12)

    def test_LiteralResponseCalculateValue6(self):
        eval_resp = EvaluateResponseExpression("4.0 - 8.0 + 3.0 * 2.0 / 5.0 ^ 2", self.optimization_problem)
        eval_resp.Initialize()
        eval_value = eval_resp.CalculateValue()
        self.assertEqual(eval_value, -3.76)

        self.assertFalse(self.optimization_problem.HasResponse(eval_resp))
        self.assertFalse(self.optimization_problem.HasResponse("4.0"))
        self.assertFalse(self.optimization_problem.HasResponse("8.0"))
        self.assertFalse(self.optimization_problem.HasResponse("3.0"))
        self.assertFalse(self.optimization_problem.HasResponse("2.0"))
        self.assertFalse(self.optimization_problem.HasResponse("5.0"))

        self.__CheckGradients([Kratos.PRESSURE, Kratos.TEMPERATURE], eval_resp, lambda : 4.0 - 8.0 + 3.0 * 2.0 / 5.0 ** 2, 1e-8, 12)

    def test_CombinedResponseCalculateValue1(self):
        eval_resp = EvaluateResponseExpression("r1", self.optimization_problem)
        eval_resp.Initialize()
        eval_value = eval_resp.CalculateValue()
        self.assertAlmostEqual(eval_value, self.r1.CalculateValue(), 12)

        self.assertTrue(self.optimization_problem.HasResponse(eval_resp))
        self.assertTrue(self.optimization_problem.HasResponse("r1"))
        self.assertEqual(eval_resp.GetName(), "r1")

        # following is the resultant response function, hence the value storage is managed by the ResponseRoutine
        self.assertFalse(ComponentDataView(self.optimization_problem.GetResponse("r1"), self.optimization_problem).GetBufferedData().HasValue("value"))

        self.__CheckGradients([Kratos.PRESSURE], eval_resp, lambda : self.r1.CalculateValue(), 1e-8, 9)

    def test_CombinedResponseCalculateValue2(self):
        eval_resp = EvaluateResponseExpression("r1 + r2 + r2", self.optimization_problem)
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

        self.__CheckGradients([Kratos.PRESSURE], eval_resp, lambda : self.r1.CalculateValue() + self.r2.CalculateValue() * 2, 1e-8, 9)

    def test_CombinedResponseCalculateValue3(self):
        eval_resp = EvaluateResponseExpression("r1 + r2 + 4.0", self.optimization_problem)
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

        self.__CheckGradients([Kratos.PRESSURE], eval_resp, lambda : self.r1.CalculateValue() + self.r2.CalculateValue() + 4.0, 1e-7, 8)

    def test_CombinedResponseCalculateValue4(self):
        eval_resp = EvaluateResponseExpression("3.5 + r1 * 2.0 + r2 / 3.0 - 4.0", self.optimization_problem)
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

        self.__CheckGradients([Kratos.PRESSURE], eval_resp, lambda : 3.5 + self.r1.CalculateValue() * 2.0 + self.r2.CalculateValue() / 3.0 - 4.0, 1e-6, 7)

    def test_CombinedResponseCalculateValue5(self):
        eval_resp = EvaluateResponseExpression("3.5 + r1 * 2.0 * r2 / 3.0 - 4.0 + r1 ^ r2", self.optimization_problem)
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

        self.__CheckGradients([Kratos.PRESSURE], eval_resp, lambda : 3.5 + self.r1.CalculateValue() * 2.0 * self.r2.CalculateValue() / 3.0 - 4.0 + self.r1.CalculateValue() ** self.r2.CalculateValue(), 1e-7, 8)

    def test_BracketResponseCalculateValue1(self):
        eval_resp = EvaluateResponseExpression("(4.0)", self.optimization_problem)
        eval_resp.Initialize()
        eval_resp.Initialize()
        eval_value = eval_resp.CalculateValue()
        self.assertEqual(eval_value, 4.0)

        with self.assertRaises(RuntimeError):
            self.__CheckGradients([Kratos.PRESSURE], eval_resp, lambda : 4.0, 1e-7, 8)

    def test_BracketResponseCalculateValue2(self):
        eval_resp = EvaluateResponseExpression("(4.0 + (6.0 / 3.0 + 3 * (2 + 8))) * 2.0", self.optimization_problem)
        eval_resp.Initialize()
        eval_value = eval_resp.CalculateValue()
        self.assertEqual(eval_value, 72.0)

        self.__CheckGradients([Kratos.PRESSURE], eval_resp, lambda : (4.0 + (6.0 / 3.0 + 3 * (2 + 8))) * 2.0, 1e-7, 8)

    def test_BracketResponseCalculateValue3(self):
        eval_resp = EvaluateResponseExpression("(4.0 + (r1 * 2) * (r2 ^ 2) + log((2+6)*(3+(4-2)/4)+6*r1) + (6.0 / 3.0 + 3 * (2 + 8))) * 2.0", self.optimization_problem)
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

        self.__CheckGradients([Kratos.PRESSURE], eval_resp, lambda : (4.0 + (self.r1.CalculateValue() * 2) * (self.r2.CalculateValue() ** 2) + log((2+6)*(3+(4-2)/4)+6*self.r1.CalculateValue()) + (6.0 / 3.0 + 3 * (2 + 8))) * 2.0, 1e-6, 7)

    def test_MultipleControlVars(self):
        eval_resp = EvaluateResponseExpression("r1 + r2 + r3 + log(r3)", self.optimization_problem)
        eval_resp.Initialize()
        self.__CheckGradients([Kratos.PRESSURE, Kratos.TEMPERATURE], eval_resp, lambda : self.r1.CalculateValue() + self.r2.CalculateValue() + self.r3.CalculateValue() + log(self.r3.CalculateValue()), 1e-8, 5)

if __name__ == "__main__":
    kratos_unittest.main()