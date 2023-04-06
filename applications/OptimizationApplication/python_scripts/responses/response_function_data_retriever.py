import KratosMultiphysics as Kratos
import KratosMultiphysics.OptimizationApplication as KratosOA
from KratosMultiphysics.OptimizationApplication.responses.response_function import ResponseFunction
from KratosMultiphysics.OptimizationApplication.utilities.optimization_info import OptimizationInfo
from KratosMultiphysics.OptimizationApplication.utilities.union_utilities import SupportedSensitivityFieldVariableTypes

class ResponseFunctionDataRetriever:
    def __init__(self, reaponse_function_name: str, optimization_info: OptimizationInfo):
        self.__response_function: ResponseFunction = optimization_info.GetOptimizationProcess(ResponseFunction, reaponse_function_name)
        self.__name = reaponse_function_name
        self.__optimization_info = optimization_info
        self.__key_prefix =  f"problem_data/response_data/{self.__name}"

    def GetPrefix(self):
        return self.__key_prefix

    def GetResponseFunction(self) -> ResponseFunction:
        return self.__response_function

    def GetScaledValue(self, step_index: int = 0, scaling: float = 1.0) -> float:
        key =f"{self.__key_prefix}/value"
        if step_index == 0 and not self.__optimization_info.HasValue(key):
            self.__optimization_info.SetValue(key, self.__response_function.CalculateValue())
        return self.__optimization_info.GetValue(key, step_index) * scaling

    def GetRelativeChange(self) -> float:
        return self.GetScaledValue() / self.GetScaledValue(1) - 1.0 if self.__optimization_info["step"] > 1 else 0.0

    def GetAbsoluteChange(self, reference_value: float) -> float:
        return self.GetScaledValue() / reference_value - 1.0 if abs(reference_value) > 1e-12 else self.GetScaledValue()

    def CalculateScaledSensitivity(self, sensitivity_variable_collective_expression_info: 'dict[SupportedSensitivityFieldVariableTypes, KratosOA.ContainerExpression.CollectiveExpressions]', scaling: float = 1.0) -> None:
        # check whether the same sensitivity is computed for each sensitivity model part w.r.t. same sensitivity variable for the same response function.
        required_sensitivity_model_part_variable_info: 'dict[SupportedSensitivityFieldVariableTypes, list[Kratos.ModelPart]]' = {}
        for sensitivity_variable, collective_expression in sensitivity_variable_collective_expression_info.items():
            is_already_computed = True
            for container_expression in collective_expression.GetContainerExpressions():
                key = f"{self.__key_prefix}/sensitivities/{container_expression.GetModelPart().FullName()}/{self.__name}_raw_{sensitivity_variable.Name()}"
                is_already_computed = is_already_computed and self.__optimization_info.HasValue(key)

            if not is_already_computed:
                # translates container expressions to model parts which is known by responses.
                required_sensitivity_model_part_variable_info[sensitivity_variable] = [container_expression.GetModelPart() for container_expression in collective_expression.GetContainerExpressions()]

        self.__response_function.CalculateSensitivity(required_sensitivity_model_part_variable_info)

        # now store the computed sensitivities from resepective containers to container expressions
        for sensitivity_variable in required_sensitivity_model_part_variable_info.keys():
            for container_expression in sensitivity_variable_collective_expression_info[sensitivity_variable].GetContainerExpressions():
                key = f"{self.__key_prefix}/sensitivities/{container_expression.GetModelPart().FullName()}/{self.__name}_raw_{sensitivity_variable.Name()}"
                container_expression.Read(sensitivity_variable)
                # this cloning is cloning of an expression pointer, hence it is not expensive.
                # no cloning of underlying sensitivity data is involved.
                self.__optimization_info.SetValue(key, container_expression.Clone())

        # now assign the sensitivities for all the collective expressions correctly
        for sensitivity_variable, collective_expression in sensitivity_variable_collective_expression_info.items():
            for container_expression in collective_expression.GetContainerExpressions():
                key = f"{self.__key_prefix}/sensitivities/{container_expression.GetModelPart().FullName()}/{self.__name}_raw_{sensitivity_variable.Name()}"
                # this copying is copying of an expression pointer, hence it is not expensive.
                # no copying of underlying sensitivity data is involved.
                container_expression.CopyFrom(self.__optimization_info.GetValue(key))
                container_expression *= scaling
