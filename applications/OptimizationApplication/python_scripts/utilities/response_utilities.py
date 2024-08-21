import re
from enum import Enum

import KratosMultiphysics as Kratos
import KratosMultiphysics.OptimizationApplication as KratosOA
from KratosMultiphysics.OptimizationApplication.responses.response_function import ResponseFunction
from KratosMultiphysics.OptimizationApplication.responses.evaluation_response_function import EvaluationResponseFunction
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

def GetClosingBracketIndex(expression: str) -> int:
    # this function assumes the starting point is without the opening bracket.
    # it will give the position of the respective closing bracket
    index = 0
    brackets_counter = 1
    while index < len(expression):
        brackets_counter += (expression[index] == "(") + (expression[index] == ")") * -1
        if brackets_counter == 0:
            return index
        index += 1

    raise RuntimeError(f"The closing bracket not found for expression ({expression}")

def GetFunction(function_name: str, optimization_problem: OptimizationProblem, *args) -> ResponseFunction:
    # import functions
    from KratosMultiphysics.OptimizationApplication.responses.log_response_function import LogResponseFunction

    functions_map: 'dict[str, type[ResponseFunction]]' = {
        "log": LogResponseFunction
    }
    if function_name in functions_map.keys():
        return functions_map[function_name](*args, optimization_problem)
    else:
        raise RuntimeError(f"Undefined \"{function_name}\" function. Followings are supported function names:\n\t" + "\n\t".join(functions_map.keys()))

def EvaluateValue(response_function: ResponseFunction, optimization_problem: OptimizationProblem) -> float:
    response_data = ComponentDataView("evaluated_responses", optimization_problem).GetBufferedData()

    if response_function.GetImplementedPhysicalKratosVariables():
        if not response_data.HasValue(f"values/{response_function.GetName()}"):
            response_data.SetValue(f"values/{response_function.GetName()}", response_function.CalculateValue())

        return response_data.GetValue(f"values/{response_function.GetName()}")
    else:
        return response_function.CalculateValue()

def EvaluateGradient(response_function: ResponseFunction, physical_variable_collective_expressions: 'dict[SupportedSensitivityFieldVariableTypes, KratosOA.CollectiveExpression]', optimization_problem: OptimizationProblem) -> 'dict[SupportedSensitivityFieldVariableTypes, KratosOA.CollectiveExpression]':
    # first get the sub_collective expressions for implemented physical kratos variables
    resp_physical_variable_collective_expressions: 'dict[SupportedSensitivityFieldVariableTypes, KratosOA.CollectiveExpression]' = {}
    for variable, collective_expression in physical_variable_collective_expressions.items():
        if variable in response_function.GetImplementedPhysicalKratosVariables():
            resp_physical_variable_collective_expressions[variable] = collective_expression.Clone()

    resp_physical_variable_collective_expressions_to_evaluate: 'dict[SupportedSensitivityFieldVariableTypes, KratosOA.CollectiveExpression]' = {}
    response_data = ComponentDataView("evaluated_responses", optimization_problem).GetUnBufferedData()

    for variable, collective_expression in resp_physical_variable_collective_expressions.items():
        # first check whether the gradients have been already evaluated. if so, take the gradients.
        if response_data.HasValue(f"gradients/d{response_function.GetName()}_d{variable.Name()}"):
            resp_physical_variable_collective_expressions[variable] = response_data.GetValue(f"gradients/d{response_function.GetName()}_d{variable.Name()}")
        else:
            # gradients have not yet evaluated. put it to the dictionary for later evaluation
            resp_physical_variable_collective_expressions_to_evaluate[variable] = collective_expression

    if len(resp_physical_variable_collective_expressions_to_evaluate) != 0:
        response_function.CalculateGradient(resp_physical_variable_collective_expressions_to_evaluate)
        for variable, collective_expression in resp_physical_variable_collective_expressions_to_evaluate.items():
            response_data.SetValue(f"gradients/d{response_function.GetName()}_d{variable.Name()}", collective_expression.Clone())
            resp_physical_variable_collective_expressions[variable] = collective_expression

    # now add zero values collective expressions to variables for which the response function does not have dependence
    for variable, collective_expression in physical_variable_collective_expressions.items():
        if variable not in response_function.GetImplementedPhysicalKratosVariables():
            resp_physical_variable_collective_expressions[variable] = collective_expression * 0.0

    return resp_physical_variable_collective_expressions

def GetResponseFunction(current_value: str, optimization_problem: OptimizationProblem) -> ResponseFunction:
    from KratosMultiphysics.OptimizationApplication.responses.literal_value_response_function import LiteralValueResponseFunction
    if current_value == "":
        return LiteralValueResponseFunction(0.0)
    else:
        try:
            # first try to get a literal response function
            return LiteralValueResponseFunction(float(current_value))
        except:
            # literal response function fails. Then the current value may hold one of the followings
            #   1. A leaf response function name
            #   2. A response expression
            #   3. A function call

            if re.match(r"(^\w+)$", current_value) is not None:
                # the match is a leaf response function name
                return optimization_problem.GetResponse(current_value)
            elif re.match(r"(^\w+\()", current_value) is not None:
                # the match is a function call
                args_starting_index = current_value.index("(")
                function_name = current_value[:args_starting_index]

                # now need to find the args correctly
                args_list: 'list[str]' = []
                index = args_starting_index + 1
                current_arg = ""
                while index < len(current_value) - 1:
                    current_char = current_value[index]

                    if current_char == "(":
                        closing_bracket_pos = GetClosingBracketIndex(current_value[index+1:])
                        current_arg += current_value[index:index+closing_bracket_pos+2]
                        index += closing_bracket_pos + 2
                        continue

                    if current_char == ",":
                        args_list.append(current_arg)
                        current_arg = ""
                    else:
                        current_arg += current_char

                    index += 1

                args_list.append(current_arg)
                return GetFunction(function_name, optimization_problem, *[__EvaluateResponseExpressionImpl(arg, optimization_problem) for arg in args_list])
            else:
                # this is a response expression.
                return __EvaluateResponseExpressionImpl(current_value, optimization_problem)

def GetValuesAndOperators(response_expression: str, optimization_problem: OptimizationProblem) -> 'tuple[list[ResponseFunction], list[str]]':
    responses: 'list[ResponseFunction]' = []
    operators: 'list[str]' = []

    index = 0
    current_word = ""
    while index < len(response_expression):
        current_char = response_expression[index]

        if current_char == "(":
            # found opening bracket
            closing_bracket_position = GetClosingBracketIndex(response_expression[index + 1:])
            if current_word != "":
                current_word = f"{current_word}{response_expression[index:index+closing_bracket_position+2]}"
            else:
                current_word = response_expression[index + 1:index+closing_bracket_position + 1]
            index += closing_bracket_position + 2
            continue

        # check for binary operator or scientific notation value
        if current_char in BinaryOperatorValues and not (re.match(r"(^[0-9.]+[e|E])$", current_word) and current_char in ["+", "-"]):
            responses.append(GetResponseFunction(current_word, optimization_problem))
            operators.append(current_char)
            current_word = ""
            index += 1
            continue

        current_word += current_char
        index += 1

    # add the last current_word
    responses.append(GetResponseFunction(current_word, optimization_problem))

    return responses, operators

def __EvaluateResponseExpressionImpl(response_expression: str, optimization_problem: OptimizationProblem) -> ResponseFunction:
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

                resultant_operator = BinaryOperatorResponseFunction(left_operand, right_operand, BinaryOperatorValuesMap[operators[operator_index]], optimization_problem)
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

def EvaluateResponseExpression(response_expression: str, optimization_problem: OptimizationProblem) -> ResponseFunction:
    evaluated_response_impl = __EvaluateResponseExpressionImpl(response_expression, optimization_problem)
    if evaluated_response_impl.GetChildResponses() and  evaluated_response_impl.GetImplementedPhysicalKratosVariables():
        # if the response has children and has some dependence on the variables, then
        # we need to use the EvaluationResponseFunction to clear the evaluation data
        # whenever CalculateValue, CalculateGradient is used.
        evaluated_response = EvaluationResponseFunction(evaluated_response_impl, optimization_problem)
        optimization_problem.AddComponent(evaluated_response)
        return evaluated_response
    else:
        # this means the following cases
        #    1. No children, but has dependence variables -> Leaf responses such as mass.
        #    2. No children, no dependence variables -> Leaf responses such as LiteralValueResponse.
        #    3. Has children, no dependence variables -> An expression with only LiteralValueResponses and functions, without any responses such as mass.
        # in the above cases, there is no need to clear the response evaluation data, hence the original evaluated
        # response is returned without the EvaluationResponseFunction wrapper.
        return evaluated_response_impl
