import KratosMultiphysics as Kratos
import KratosMultiphysics.OptimizationApplication as KratosOA
from KratosMultiphysics.OptimizationApplication.utilities.optimization_info import OptimizationInfo
from KratosMultiphysics.OptimizationApplication.responses.response_function import ResponseFunction
from KratosMultiphysics.OptimizationApplication.responses.response_utilities import CalculateRelativeChange
from KratosMultiphysics.OptimizationApplication.responses.response_utilities import CalculateResponseSensitivity
from KratosMultiphysics.OptimizationApplication.responses.response_utilities import GetResponseValue
from KratosMultiphysics.OptimizationApplication.utilities.union_utilities import SupportedSensitivityFieldVariableTypes

class StandardizedConstraint:
    def __init__(self, parameters: Kratos.Parameters, optimization_info: OptimizationInfo):
        default_parameters = Kratos.Parameters("""{
            "response_name": "",
            "type"         : "",
            "scaling"      : 1.0,
            "ref_type"     : "",
            "ref_value"    : 0.0
        }""")
        parameters.ValidateAndAssignDefaults(default_parameters)

        self.__scaling = parameters["scaling"].GetDouble()

        self.__ref_type = parameters["ref_type"].GetString()
        if self.__ref_type == "initial_value":
            self.__reference_value = None
        elif self.__ref_type == "specified_value":
            self.__reference_value = parameters["ref_value"].GetDouble()
        else:
            raise RuntimeError(f"Provided \"reference_type\" = {self.__ref_type} is not supported for constraint response functions. Followings are supported options: \n\tinitial_value\n\tspecified_value")

        self.__constraint_type = parameters["type"].GetString()
        if self.__constraint_type in ["<", "<=", "="]:
            self.__standardization_value = 1.0
        elif self.__constraint_type in [">", ">="]:
            self.__standardization_value = -1.0
        else:
            raise RuntimeError(f"Provided \"type\" = {self.__constraint_type} is not supported in constraint response functions. Followings are supported options: \n\t=\n\t<=\n\t<\n\t>=\n\t>")

        self.__name = parameters["response_name"].GetString()
        self.__response_function: ResponseFunction = optimization_info.GetOptimizationProcess(ResponseFunction, self.__name)

        self.__initial_response_value = None
        self.__optimization_info = optimization_info
        self._key_prefix = f"problem_data/response_data/{self.GetName()}"

    def GetName(self) -> str:
        return self.__name

    def GetResponseType(self) -> str:
        return self.__constraint_type

    def GetResponseFunction(self) -> ResponseFunction:
        return self.__response_function

    def GetReferenceValue(self):
        if self.__reference_value is None:
            self.__reference_value = self.GetValue()

        return self.__reference_value

    def GetValue(self, step_index: int = 0) -> float:
        return GetResponseValue(self.GetName(), self.__optimization_info, step_index)

    def GetStandardizedValue(self, step_index: int = 0) -> float:
        return self.__standardization_value * (self.GetValue(step_index) * self.__scaling - self.GetReferenceValue())

    def IsActive(self, step_index: int = 0) -> bool:
        return (self.__constraint_type == "=") or self.GetStandardizedValue(step_index) >= 0.0

    def GetViolationRatio(self, step_index: int = 0) -> float:
        if self.IsActive():
            return CalculateRelativeChange(self.__standardization_value * self.GetValue(step_index) * self.__scaling, self.GetReferenceValue())
        else:
            return 0.0

    def CalculateSensitivity(self, sensitivity_variable_collective_expression_info: 'dict[SupportedSensitivityFieldVariableTypes, KratosOA.ContainerExpression.CollectiveExpressions]'):
        CalculateResponseSensitivity(self.GetName(), self.__optimization_info, sensitivity_variable_collective_expression_info)

    def CalculateStandardizedSensitivity(self, sensitivity_variable_collective_expression_info: 'dict[SupportedSensitivityFieldVariableTypes, KratosOA.ContainerExpression.CollectiveExpressions]'):
        self.CalculateSensitivity(sensitivity_variable_collective_expression_info)
        for value in sensitivity_variable_collective_expression_info.values():
            value *= (self.__standardization_value * self.__scaling)

    def UpdateOptimizationInfo(self) -> None:
        absolute_change = CalculateRelativeChange(self.GetValue(), self.__GetInitialResponseValue())

        if self.__optimization_info["step"] > 1:
            relative_change = CalculateRelativeChange(self.GetValue(), self.GetValue(1))
        else:
            relative_change = 0.0

        self.__optimization_info.SetValue(f"{self._key_prefix}/relative_change", relative_change)
        self.__optimization_info.SetValue(f"{self._key_prefix}/absolute_change", absolute_change)
        self.__optimization_info.SetValue(f"{self._key_prefix}/reference_value", self.GetReferenceValue())
        self.__optimization_info.SetValue(f"{self._key_prefix}/is_active", self.IsActive())
        self.__optimization_info.SetValue(f"{self._key_prefix}/type", self.GetResponseType())
        self.__optimization_info.SetValue(f"{self._key_prefix}/violation_ratio", abs(self.GetViolationRatio()))

    def GetResponseInfo(self) -> str:
        if not (self.__optimization_info.HasValue(f"{self._key_prefix}/relative_change")):
            self.UpdateOptimizationInfo()

        relative_change = self.__optimization_info.GetValue(f"{self._key_prefix}/relative_change")
        absolute_change = self.__optimization_info.GetValue(f"{self._key_prefix}/absolute_change")
        violation = self.__optimization_info.GetValue(f"{self._key_prefix}/violation_ratio")

        msg = "\tConstraint info:"
        msg += f"\n\t\t name          : {self.GetName()}"
        msg += f"\n\t\t value         : {self.GetValue()}"
        msg += f"\n\t\t type          : {self.GetResponseType()}"
        msg += f"\n\t\t ref_value     : {self.GetReferenceValue()}"
        msg += f"\n\t\t is_active     : {self.IsActive()}"
        msg += f"\n\t\t rel_change [%]: {relative_change * 100.0}"
        msg += f"\n\t\t abs_change [%]: {absolute_change * 100.0}"
        msg += f"\n\t\t violation  [%]: {violation * 100.0}"

        return msg

    def __GetInitialResponseValue(self):
        if self.__initial_response_value is None:
            self.__initial_response_value = self.GetValue()

        return self.__initial_response_value