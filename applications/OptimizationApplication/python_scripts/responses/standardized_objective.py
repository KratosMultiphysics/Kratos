import KratosMultiphysics as Kratos
import KratosMultiphysics.OptimizationApplication as KratosOA
from KratosMultiphysics.OptimizationApplication.utilities.optimization_info import OptimizationInfo
from KratosMultiphysics.OptimizationApplication.responses.response_function import ResponseFunction
from KratosMultiphysics.OptimizationApplication.responses.response_utilities import CalculateRelativeChange
from KratosMultiphysics.OptimizationApplication.responses.response_utilities import CalculateResponseSensitivity
from KratosMultiphysics.OptimizationApplication.responses.response_utilities import GetResponseValue
from KratosMultiphysics.OptimizationApplication.utilities.union_utilities import SupportedSensitivityFieldVariableTypes

class StandardizedObjective:
    def __init__(self, parameters: Kratos.Parameters, optimization_info: OptimizationInfo):
        default_parameters = Kratos.Parameters("""{
            "response_name": "",
            "type"         : "",
            "scaling"      : 1.0
        }""")
        parameters.ValidateAndAssignDefaults(default_parameters)

        self.__scaling = parameters["scaling"].GetDouble()

        self.__objective_type = parameters["type"].GetString()
        if self.__objective_type == "minimization":
            self.__scaling *= 1.0
        elif self.__objective_type == "maximization":
            self.__scaling *= -1.0
        else:
            raise RuntimeError(f"Requesting unsupported type {self.__objective_type} for objective response function. Supported types are: \n\tminimization\n\tmaximization")

        self.__name = parameters["response_name"].GetString()
        self.__response_function: ResponseFunction = optimization_info.GetOptimizationProcess(ResponseFunction, self.__name)

        self.__initial_response_value = None
        self.__optimization_info = optimization_info
        self._key_prefix = f"problem_data/response_data/{self.GetName()}"

    def GetName(self) -> str:
        return self.__name

    def GetResponseFunction(self) -> ResponseFunction:
        return self.__response_function

    def GetResponseType(self) -> str:
        return self.__objective_type

    def GetValue(self, step_index: int = 0) -> float:
        return GetResponseValue(self.GetName(), self.__optimization_info, step_index)

    def GetStandardizedValue(self, step_index: int = 0) -> float:
        return self.GetValue(step_index) * self.__scaling

    def CalculateSensitivity(self, sensitivity_variable_collective_expression_info: 'dict[SupportedSensitivityFieldVariableTypes, KratosOA.ContainerExpression.CollectiveExpressions]'):
        CalculateResponseSensitivity(self.GetName(), self.__optimization_info, sensitivity_variable_collective_expression_info)

    def CalculateStandardizedSensitivity(self, sensitivity_variable_collective_expression_info: 'dict[SupportedSensitivityFieldVariableTypes, KratosOA.ContainerExpression.CollectiveExpressions]'):
        self.CalculateSensitivity(sensitivity_variable_collective_expression_info)
        for value in sensitivity_variable_collective_expression_info.values():
            value *= self.__scaling

    def UpdateOptimizationInfo(self) -> None:
        absolute_change = CalculateRelativeChange(self.GetValue(), self.__GetInitialResponseValue())
        if self.__optimization_info["step"] > 1:
            relative_change = CalculateRelativeChange(self.GetValue(), self.GetValue(1))
        else:
            relative_change = 0.0
        self.__optimization_info.SetValue(f"{self._key_prefix}/type", self.GetResponseType())
        self.__optimization_info.SetValue(f"{self._key_prefix}/relative_change", relative_change)
        self.__optimization_info.SetValue(f"{self._key_prefix}/absolute_change", absolute_change)

    def GetResponseInfo(self) -> str:
        if not (self.__optimization_info.HasValue(f"{self._key_prefix}/relative_change")):
            self.UpdateOptimizationInfo()

        relative_change = self.__optimization_info.GetValue(f"{self._key_prefix}/relative_change")
        absolute_change = self.__optimization_info.GetValue(f"{self._key_prefix}/absolute_change")
        msg = "\tObjective info:"
        msg += f"\n\t\t name          : {self.GetName()}"
        msg += f"\n\t\t type          : {self.GetResponseType()}"
        msg += f"\n\t\t value         : {self.GetValue()}"
        msg += f"\n\t\t rel_change [%]: {relative_change * 100.0}"
        msg += f"\n\t\t abs_change [%]: {absolute_change * 100.0}"
        return msg

    def __GetInitialResponseValue(self):
        if self.__initial_response_value is None:
            self.__initial_response_value = self.GetValue()

        return self.__initial_response_value