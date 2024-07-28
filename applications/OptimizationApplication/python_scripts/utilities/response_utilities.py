from enum import Enum

import KratosMultiphysics as Kratos
import KratosMultiphysics.OptimizationApplication as KratosOA
from KratosMultiphysics.OptimizationApplication.responses.response_function import ResponseFunction
from KratosMultiphysics.OptimizationApplication.responses.response_function import SupportedSensitivityFieldVariableTypes
from KratosMultiphysics.OptimizationApplication.utilities.union_utilities import SupportedSensitivityFieldVariableTypes
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem import OptimizationProblem
from KratosMultiphysics.OptimizationApplication.utilities.component_data_view import ComponentDataView

class BinaryOperator(str, Enum):
    POWER    = "^"
    DIVIDE   = "/"
    MULTIPLY = "*"
    ADD      = "+"
    SUBTRACT = "-"

BinaryOperatorValues = [operator.value for operator in BinaryOperator]

BinaryOperatorValuesMap = dict([(operator.value, operator) for operator in BinaryOperator])

def GetFunctionsMap() -> 'dict[str, ResponseFunction]':
    # import functions
    from KratosMultiphysics.OptimizationApplication.responses.log_response_function import LogResponseFunction
    return {
        "log": LogResponseFunction
    }

def EvaluateValue(response_function: ResponseFunction, optimization_problem: OptimizationProblem) -> float:
    if optimization_problem.HasResponse(response_function):
        response_data = ComponentDataView(response_function, optimization_problem)
        if response_data.HasDataBuffer():
            response_data_buffer = response_data.GetBufferedData()
            if not response_data_buffer.HasValue("value"):
                response_data_buffer.SetValue("value", response_function.CalculateValue())
            return response_data_buffer.GetValue("value")
        else:
            return response_function.CalculateValue()
    else:
        return response_function.CalculateValue()

def EvaluateGradient(response_function: ResponseFunction, physical_variable_collective_expressions: 'dict[SupportedSensitivityFieldVariableTypes, KratosOA.CollectiveExpression]', optimization_problem: OptimizationProblem) -> 'dict[SupportedSensitivityFieldVariableTypes, KratosOA.CollectiveExpression]':
    # first get the sub_collective expressions for implemented physical kratos variables
    resp_physical_variable_collective_expressions: 'dict[SupportedSensitivityFieldVariableTypes, KratosOA.CollectiveExpression]' = {}
    for variable, collective_expression in physical_variable_collective_expressions.items():
        if variable in response_function.GetImplementedPhysicalKratosVariables():
            resp_physical_variable_collective_expressions[variable] = collective_expression.Clone()

    resp_physical_variable_collective_expressions_to_evaluate: 'dict[SupportedSensitivityFieldVariableTypes, KratosOA.CollectiveExpression]' = {}
    response_data = ComponentDataView(response_function, optimization_problem)
    if response_data.HasDataBuffer():
        response_data_buffer = response_data.GetBufferedData()
        for variable, collective_expression in resp_physical_variable_collective_expressions.items():
            # first check whether the gradients have been already evaluated. if so, take the gradients.
            if response_data_buffer.HasValue(f"d{response_function.GetName()}_d{variable.Name()}"):
                resp_physical_variable_collective_expressions[variable] = response_data_buffer.GetValue(f"d{response_function.GetName()}_d{variable.Name()}")
            else:
                # gradients have not yet evaluated. put it to the dictionary for later evaluation
                resp_physical_variable_collective_expressions_to_evaluate[variable] = collective_expression
    else:
        for variable, collective_expression in resp_physical_variable_collective_expressions.items():
            # gradients have not yet evaluated. put it to the dictionary for later evaluation
            resp_physical_variable_collective_expressions_to_evaluate[variable] = collective_expression

    if len(resp_physical_variable_collective_expressions_to_evaluate) != 0:
        response_function.CalculateGradient(resp_physical_variable_collective_expressions_to_evaluate)
        for variable, collective_expression in resp_physical_variable_collective_expressions_to_evaluate.items():
            if response_data.HasDataBuffer():
                response_data.GetBufferedData().SetValue(f"d{response_function.GetName()}_d{variable.Name()}", collective_expression.Clone())
            resp_physical_variable_collective_expressions[variable] = collective_expression

    # now add zero values collective expressions to variables for which the response function does not have dependence
    for variable, collective_expression in physical_variable_collective_expressions.items():
        if variable not in response_function.GetImplementedPhysicalKratosVariables():
            resp_physical_variable_collective_expressions[variable] = collective_expression * 0.0

    return resp_physical_variable_collective_expressions

def GetResponseFunction(current_value: str, optimization_problem: OptimizationProblem) -> ResponseFunction:
    from KratosMultiphysics.OptimizationApplication.responses.literal_value_response_function import LiteralValueResponseFunction
    list_of_response_names: 'list[str]' = [response.GetName() for response in optimization_problem.GetListOfResponses()]
    if current_value == "":
        return LiteralValueResponseFunction(0.0)
    else:
        try:
            return LiteralValueResponseFunction(float(current_value))
        except:
            # first check whether the value exists in list of responses
            if current_value in list_of_response_names:
                return optimization_problem.GetResponse(current_value)
            else:
                raise RuntimeError(f"The response named \"{current_value}\" not defined.")

def GetValuesAndOperators(response_expression: str, optimization_problem: OptimizationProblem) -> 'tuple[list[ResponseFunction], list[str]]':
    responses: 'list[ResponseFunction]' = []
    operators: 'list[str]' = []

    index = 0
    current_word = ""
    while index < len(response_expression):
        current_char = response_expression[index]

        if current_char in BinaryOperatorValues:
            responses.append(GetResponseFunction(current_word, optimization_problem))
            operators.append(current_char)
            current_word = ""
        else:
            current_word += current_char

        index += 1

    # add the last current_word
    responses.append(GetResponseFunction(current_word, optimization_problem))

    return responses, operators

def EvaluateResponseExpression(model: Kratos.Model, response_expression: str, optimization_problem: OptimizationProblem) -> ResponseFunction:
    from KratosMultiphysics.OptimizationApplication.responses.binary_operator_response_function import BinaryOperatorResponseFunction

    response_expression = response_expression.replace(" ", "")
    responses, operators = GetValuesAndOperators(response_expression, optimization_problem)

    def __evaluate_operator(list_of_operators: 'list[str]') -> None:
        operator_index = 0
        while operator_index < len(operators):
            if operators[operator_index] in list_of_operators:
                left_operand  = responses[operator_index]
                right_operand = responses[operator_index + 1]

                if not optimization_problem.HasResponse(left_operand) and left_operand.GetImplementedPhysicalKratosVariables() != []:
                    # the left operand is not in the optimization, and it will not be accessible since it will
                    # be within the BinaryOperatorResponseFunction. Hence adding it to optimization problem
                    # so that we can put the intermediate values and gradients to the optimization problem.
                    optimization_problem.AddComponent(left_operand)

                if not optimization_problem.HasResponse(right_operand) and right_operand.GetImplementedPhysicalKratosVariables() != []:
                    # the right operand is not in the optimization, and it will not be accessible since it will
                    # be within the BinaryOperatorResponseFunction. Hence adding it to optimization problem
                    # so that we can put the intermediate values and gradients to the optimization problem.
                    optimization_problem.AddComponent(right_operand)

                # now check whether the left operand is there in the optimization problem
                if optimization_problem.HasResponse(left_operand):
                    left_operand_component_data_view = ComponentDataView(left_operand, optimization_problem)
                    # if it is in the optimization problem, and no data buffer has been defined, then define one.
                    if not left_operand_component_data_view.HasDataBuffer():
                        left_operand_component_data_view.SetDataBuffer(1)

                # now check whether the right operand is there in the optimization problem
                if optimization_problem.HasResponse(right_operand):
                    right_operand_component_data_view = ComponentDataView(right_operand, optimization_problem)
                    # if it is in the optimization problem, and no data buffer has been defined, then define one.
                    if not right_operand_component_data_view.HasDataBuffer():
                        right_operand_component_data_view.SetDataBuffer(1)

                resultant_operator = BinaryOperatorResponseFunction(model, left_operand, right_operand, BinaryOperatorValuesMap[operators[operator_index]], optimization_problem)
                if not optimization_problem.HasResponse(resultant_operator.GetName()):
                    # There is no already existing response with the same operation. Then use the new one.
                    responses[operator_index] = resultant_operator
                else:
                    # There is an already existing response with the same operation. Then take it.
                    responses[operator_index] = optimization_problem.GetResponse(resultant_operator.GetName())

                del responses[operator_index + 1]
                del operators[operator_index]
            else:
                operator_index += 1

    # according to BODMAS
    __evaluate_operator(["^"])
    __evaluate_operator(["*", "/"])
    __evaluate_operator(["+", "-"])

    return responses[0]