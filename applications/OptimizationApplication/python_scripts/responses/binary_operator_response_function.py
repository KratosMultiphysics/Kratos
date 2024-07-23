from typing import Optional
from enum import Enum
from math import log

import KratosMultiphysics as Kratos
import KratosMultiphysics.OptimizationApplication as KratosOA
from KratosMultiphysics.OptimizationApplication.responses.response_function import ResponseFunction
from KratosMultiphysics.OptimizationApplication.responses.response_function import SupportedSensitivityFieldVariableTypes
from KratosMultiphysics.OptimizationApplication.utilities.union_utilities import SupportedSensitivityFieldVariableTypes
from KratosMultiphysics.OptimizationApplication.utilities.model_part_utilities import ModelPartOperation
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem import OptimizationProblem
from KratosMultiphysics.OptimizationApplication.utilities.component_data_view import ComponentDataView
from KratosMultiphysics.OptimizationApplication.utilities.response_utilities import EvaluateGradient

class BinaryOperatorResponseFunction(ResponseFunction):
    class BinaryOperator(str, Enum):
        ADD      = "+"
        SUBTRACT = "-"
        MULTIPLY = "*"
        DIVIDE   = "/"
        POWER    = "^"

    def __init__(self, model: Kratos.Model, response_function_1: ResponseFunction, response_function_2: ResponseFunction, binary_operator: BinaryOperator, optimization_problem: OptimizationProblem):
        super().__init__(f"({response_function_1.GetName()}{binary_operator.value}{response_function_2.GetName()})")

        if response_function_1 not in optimization_problem.GetListOfResponses():
            optimization_problem.AddComponent(response_function_1)

        if response_function_2 not in optimization_problem.GetListOfResponses():
            optimization_problem.AddComponent(response_function_2)

        self.model = model
        self.optimization_problem = optimization_problem
        self.response_function_1 = response_function_1
        self.response_function_2 = response_function_2
        self.binary_operator = binary_operator
        self.model_part: Optional[Kratos.ModelPart] = None

    def GetImplementedPhysicalKratosVariables(self) -> 'list[SupportedSensitivityFieldVariableTypes]':
        vars_list = self.response_function_1.GetImplementedPhysicalKratosVariables()
        vars_list.extend(self.response_function_2.GetImplementedPhysicalKratosVariables())
        return vars_list

    def Initialize(self) -> None:
        self.response_function_1.Initialize()
        self.response_function_2.Initialize()
        self.model_part = ModelPartOperation(self.model, ModelPartOperation.OperationType.UNION, f"response_{self.GetName()}", [self.response_function_1.GetInfluencingModelPart().FullName(), self.response_function_2.GetInfluencingModelPart().FullName()], False).GetModelPart()

    def Check(self) -> None:
        self.response_function_1.Check()
        self.response_function_2.Check()

    def Finalize(self) -> None:
        self.response_function_1.Finalize()
        self.response_function_2.Finalize()

    def GetInfluencingModelPart(self) -> Kratos.ModelPart:
        return self.model_part

    def CalculateValue(self) -> float:
        # check whether the response function 1 is already evaluated for current step
        resp_1_buffered_data = ComponentDataView(self.response_function_1, self.optimization_problem).GetBufferedData()
        if resp_1_buffered_data.HasValue("value"):
            # if it is evaluated, then take the value
            resp_1_value = resp_1_buffered_data.GetValue("value")
        else:
            # if not, compute the value
            resp_1_value = self.response_function_1.CalculateValue()
            resp_1_buffered_data.SetValue("value", resp_1_value)

        # check whether the response function 2 is already evaluated for current step
        resp_2_buffered_data = ComponentDataView(self.response_function_2, self.optimization_problem).GetBufferedData()
        if resp_2_buffered_data.HasValue("value"):
            # if it is evaluated, then take the value
            resp_2_value = resp_2_buffered_data.GetValue("value")
        else:
            # if not, compute the value
            resp_2_value = self.response_function_1.CalculateValue()
            resp_2_buffered_data.SetValue("value", resp_2_value)

        # now do the binary arithmetics.
        if self.binary_operator == BinaryOperatorResponseFunction.BinaryOperator.ADD:
            return resp_1_value + resp_2_value
        elif self.binary_operator == BinaryOperatorResponseFunction.BinaryOperator.SUBTRACT:
            return resp_1_value - resp_2_value
        elif self.binary_operator == BinaryOperatorResponseFunction.BinaryOperator.MULTIPLY:
            return resp_1_value * resp_2_value
        elif self.binary_operator == BinaryOperatorResponseFunction.BinaryOperator.DIVIDE:
            return resp_1_value / resp_2_value
        elif self.binary_operator == BinaryOperatorResponseFunction.BinaryOperator.POWER:
            return resp_1_value ** resp_2_value

    def CalculateGradient(self, physical_variable_collective_expressions: 'dict[SupportedSensitivityFieldVariableTypes, KratosOA.CollectiveExpression]') -> None:
        resp_1_buffered_data = ComponentDataView(self.response_function_1, self.optimization_problem).GetBufferedData()
        resp_2_buffered_data = ComponentDataView(self.response_function_2, self.optimization_problem).GetBufferedData()

        resp_1_physical_variable_collective_expressions = EvaluateGradient(self.response_function_1, resp_1_buffered_data, physical_variable_collective_expressions)
        resp_2_physical_variable_collective_expressions = EvaluateGradient(self.response_function_2, resp_2_buffered_data, physical_variable_collective_expressions)
        v1: float = resp_1_buffered_data.GetValue("value")
        v2: float = resp_2_buffered_data.GetValue("value")

        for variable, collective_expression in physical_variable_collective_expressions.items():
            for result, g1, g2 in zip(collective_expression.GetContainerExpressions(), resp_1_physical_variable_collective_expressions[variable].GetContainerExpressions(), resp_2_physical_variable_collective_expressions[variable].GetContainerExpressions()):
                if self.binary_operator == BinaryOperatorResponseFunction.BinaryOperator.ADD:
                    result.SetExpression((g1 + g2).GetExpression())
                elif self.binary_operator == BinaryOperatorResponseFunction.BinaryOperator.SUBTRACT:
                    result.SetExpression((g1 - g2).GetExpression())
                elif self.binary_operator == BinaryOperatorResponseFunction.BinaryOperator.MULTIPLY:
                    result.SetExpression((g1 * v2 + g2 * v1).GetExpression())
                elif self.binary_operator == BinaryOperatorResponseFunction.BinaryOperator.DIVIDE:
                    result.SetExpression((g1 / v2 - g2 * (v1 / v2 ** 2)).GetExpression())
                elif self.binary_operator == BinaryOperatorResponseFunction.BinaryOperator.POWER:
                    result.SetExpression(((g1 * (v2 / v1) + g2 * log(v1)) * (v1 ** v2)).GetExpression())

    def __str__(self) -> str:
        return f"Response [type = {self.__class__.__name__}, name = {self.GetName()}, model part name = {self.model_part.FullName()}]"