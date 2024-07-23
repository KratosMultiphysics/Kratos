import KratosMultiphysics.OptimizationApplication as KratosOA
from KratosMultiphysics.OptimizationApplication.responses.response_function import ResponseFunction
from KratosMultiphysics.OptimizationApplication.responses.response_function import SupportedSensitivityFieldVariableTypes
from KratosMultiphysics.OptimizationApplication.utilities.union_utilities import SupportedSensitivityFieldVariableTypes
from KratosMultiphysics.OptimizationApplication.utilities.buffered_dict import BufferedDict

def EvaluateGradient(response_function: ResponseFunction, response_buffered_data: BufferedDict, physical_variable_collective_expressions: 'dict[SupportedSensitivityFieldVariableTypes, KratosOA.CollectiveExpression]') -> 'dict[SupportedSensitivityFieldVariableTypes, KratosOA.CollectiveExpression]':
    # first get the sub_collective expressions for implemented physical kratos variables
    resp_physical_variable_collective_expressions: 'dict[SupportedSensitivityFieldVariableTypes, KratosOA.CollectiveExpression]' = {}
    for variable, collective_expression in physical_variable_collective_expressions.items():
        if variable in response_function.GetImplementedPhysicalKratosVariables():
            resp_physical_variable_collective_expressions[variable] = collective_expression.Clone()

    resp_physical_variable_collective_expressions_to_evaluate: 'dict[SupportedSensitivityFieldVariableTypes, KratosOA.CollectiveExpression]' = {}
    for variable, collective_expression in resp_physical_variable_collective_expressions.items():
        # first check whether the gradients have been already evaluated. if so, take the gradients.
        if response_buffered_data.HasValue(f"d{response_function.GetName()}_d{variable.Name()}"):
            resp_physical_variable_collective_expressions[variable] = response_buffered_data.GetValue(f"d{response_function.GetName()}_d{variable.Name()}")
        else:
            # gradients have not yet evaluated. put it to the dictionary for later evaluation
            resp_physical_variable_collective_expressions_to_evaluate[variable] = collective_expression

    if len(resp_physical_variable_collective_expressions_to_evaluate) != 0:
        response_function.CalculateGradient(resp_physical_variable_collective_expressions_to_evaluate)
        for variable, collective_expression in resp_physical_variable_collective_expressions_to_evaluate.items():
            response_buffered_data.SetValue(f"d{response_function.GetName()}_d{variable.Name()}", collective_expression.Clone())
            resp_physical_variable_collective_expressions[variable] = collective_expression

    # now add zero values collective expressions to variables for which the response function does not have dependence
    for variable, collective_expression in physical_variable_collective_expressions.items():
        if variable not in response_function.GetImplementedPhysicalKratosVariables():
            resp_physical_variable_collective_expressions[variable] = collective_expression * 0.0

    return resp_physical_variable_collective_expressions