import KratosMultiphysics as Kratos
import KratosMultiphysics.OptimizationApplication as KratosOA
from KratosMultiphysics.OptimizationApplication.responses.response_function import ResponseFunction
from KratosMultiphysics.OptimizationApplication.utilities.optimization_info import OptimizationInfo
from KratosMultiphysics.OptimizationApplication.utilities.union_utilities import SupportedSensitivityFieldVariableTypes

def CalculateRelativeChange(current_value: float, previous_value: float, tolerance = 1e-12) -> float:
    if abs(previous_value) < tolerance:
        return current_value
    else:
        return (current_value / previous_value - 1.0)

def GetResponseValue(response_name: str, optimization_info: OptimizationInfo, step_index = 0) -> float:
    response_function: ResponseFunction = optimization_info.GetOptimizationProcess(ResponseFunction, response_name)
    key = f"problem_data/response_data/{response_name}/value"
    if step_index == 0 and not optimization_info.HasValue(key):
        optimization_info.SetValue(key, response_function.CalculateValue())
    return optimization_info.GetValue(key, step_index)

def CalculateResponseSensitivity(response_name: str, optimization_info: OptimizationInfo, sensitivity_variable_collective_expression_info: 'dict[SupportedSensitivityFieldVariableTypes, KratosOA.ContainerExpression.CollectiveExpressions]') -> None:
    response_function: ResponseFunction = optimization_info.GetOptimizationProcess(ResponseFunction, response_name)
    key_prefix = f"problem_data/response_data/{response_name}"

    # check whether the same sensitivity is computed for each sensitivity model part w.r.t. same sensitivity variable for the same response function.
    required_sensitivity_model_part_variable_info: 'dict[SupportedSensitivityFieldVariableTypes, list[Kratos.ModelPart]]' = {}
    for sensitivity_variable, collective_expression in sensitivity_variable_collective_expression_info.items():
        is_already_computed = True
        for container_expression in collective_expression.GetContainerExpressions():
            key = f"{key_prefix}/sensitivities/{container_expression.GetModelPart().FullName()}/{response_name}_{sensitivity_variable.Name()}_raw"
            is_already_computed = is_already_computed and optimization_info.HasValue(key)

        if not is_already_computed:
            # translates container expressions to model parts which is known by responses.
            required_sensitivity_model_part_variable_info[sensitivity_variable] = [container_expression.GetModelPart() for container_expression in collective_expression.GetContainerExpressions()]

    response_function.CalculateSensitivity(required_sensitivity_model_part_variable_info)

    # now store the computed sensitivities from resepective containers to container expressions
    for sensitivity_variable in required_sensitivity_model_part_variable_info.keys():
        for container_expression in sensitivity_variable_collective_expression_info[sensitivity_variable].GetContainerExpressions():
            key = f"{key_prefix}/sensitivities/{container_expression.GetModelPart().FullName()}/{response_name}_{sensitivity_variable.Name()}_raw"
            container_expression.Read(sensitivity_variable)
            # this cloning is cloning of an expression pointer, hence it is not expensive.
            # no cloning of underlying sensitivity data is involved.
            optimization_info.SetValue(key, container_expression.Clone())

    # now assign the sensitivities for all the collective expressions correctly
    for sensitivity_variable, collective_expression in sensitivity_variable_collective_expression_info.items():
        for container_expression in collective_expression.GetContainerExpressions():
            key = f"{key_prefix}/sensitivities/{container_expression.GetModelPart().FullName()}/{response_name}_{sensitivity_variable.Name()}_raw"
            # this copying is copying of an expression pointer, hence it is not expensive.
            # no copying of underlying sensitivity data is involved.
            container_expression.CopyFrom(optimization_info.GetValue(key))