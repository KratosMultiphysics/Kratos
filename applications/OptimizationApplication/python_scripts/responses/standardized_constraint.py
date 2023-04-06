import KratosMultiphysics as Kratos
import KratosMultiphysics.OptimizationApplication as KratosOA
from KratosMultiphysics.OptimizationApplication.utilities.optimization_info import OptimizationInfo
from KratosMultiphysics.OptimizationApplication.responses.response_function_data_retriever import ResponseFunctionDataRetriever
from KratosMultiphysics.OptimizationApplication.utilities.union_utilities import SupportedSensitivityFieldVariableTypes

class StandardizedConstraint:
    def __init__(self, parameters: Kratos.Parameters, optimization_info: OptimizationInfo):
        default_parameters = Kratos.Parameters("""{
            "response_name": "",
            "type"         : "",
            "ref_value"    : 0.0
        }""")

        if parameters.Has("ref_value") and parameters.IsString("ref_value"):
            default_parameters["ref_value"].SetString("")

        parameters.ValidateAndAssignDefaults(default_parameters)

        if parameters["ref_value"].IsDouble():
            self.__ref_type = "specified_value"
            self.__reference_value = parameters["ref_value"].GetDouble()
        elif parameters["ref_value"].IsString() and parameters["ref_value"].GetString() == "initial_value":
            self.__ref_type = "initial_value"
            self.__reference_value = None
        else:
            raise RuntimeError(f"Provided \"reference_type\" = {self.__ref_type} is not supported for constraint response functions. Followings are supported options: \n\tinitial_value\n\tfloat value")

        self.__constraint_type = parameters["type"].GetString()
        if self.__constraint_type in ["<", "<=", "="]:
            self.__standardization_value = 1.0
        elif self.__constraint_type in [">", ">="]:
            self.__standardization_value = -1.0
        else:
            raise RuntimeError(f"Provided \"type\" = {self.__constraint_type} is not supported in constraint response functions. Followings are supported options: \n\t=\n\t<=\n\t<\n\t>=\n\t>")

        self.__name = parameters["response_name"].GetString()
        self.__response_function_data_retriever = ResponseFunctionDataRetriever(self.GetName(), optimization_info)
        self.__optimization_info = optimization_info

    def GetName(self) -> str:
        return self.__name

    def GetResponseType(self) -> str:
        return self.__constraint_type

    def GetReferenceValue(self):
        if self.__reference_value is None:
            self.__reference_value = self.__response_function_data_retriever.GetScaledValue()
        return self.__reference_value

    def GetStandardizedValue(self, step_index: int = 0) -> float:
        return self.__response_function_data_retriever.GetScaledValue(step_index, self.__standardization_value) - self.__standardization_value * self.GetReferenceValue()

    def IsActive(self, step_index: int = 0) -> bool:
        return (self.__constraint_type == "=") or self.GetStandardizedValue(step_index) >= 0.0

    def GetViolationRatio(self, step_index: int = 0) -> float:
        return (self.GetStandardizedValue(step_index) / self.GetReferenceValue() if abs(self.GetReferenceValue()) > 1e-12 else self.GetStandardizedValue(step_index)) * self.IsActive()

    def CalculateStandardizedSensitivity(self, sensitivity_variable_collective_expression_info: 'dict[SupportedSensitivityFieldVariableTypes, KratosOA.ContainerExpression.CollectiveExpressions]'):
        self.__response_function_data_retriever.CalculateScaledSensitivity(sensitivity_variable_collective_expression_info, self.__standardization_value)

    def UpdateOptimizationInfo(self) -> None:
        self.__optimization_info.SetValue(f"{self.__response_function_data_retriever.GetPrefix()}/relative_change", self.__response_function_data_retriever.GetRelativeChange())
        self.__optimization_info.SetValue(f"{self.__response_function_data_retriever.GetPrefix()}/absolute_change", self.__response_function_data_retriever.GetAbsoluteChange(self.GetReferenceValue()))
        self.__optimization_info.SetValue(f"{self.__response_function_data_retriever.GetPrefix()}/reference_value", self.GetReferenceValue())
        self.__optimization_info.SetValue(f"{self.__response_function_data_retriever.GetPrefix()}/is_active", self.IsActive())
        self.__optimization_info.SetValue(f"{self.__response_function_data_retriever.GetPrefix()}/type", self.GetResponseType())
        self.__optimization_info.SetValue(f"{self.__response_function_data_retriever.GetPrefix()}/violation_ratio", abs(self.GetViolationRatio()))

    def GetResponseInfo(self) -> str:
        msg = "\tConstraint info:"
        msg += f"\n\t\t name          : {self.GetName()}"
        msg += f"\n\t\t value         : {self.GetValue()}"
        msg += f"\n\t\t type          : {self.GetResponseType()}"
        msg += f"\n\t\t ref_value     : {self.GetReferenceValue()}"
        msg += f"\n\t\t is_active     : {self.IsActive()}"
        msg += f"\n\t\t rel_change [%]: {self.__response_function_data_retriever.GetRelativeChange() * 100.0}"
        msg += f"\n\t\t abs_change [%]: {self.__response_function_data_retriever.GetAbsoluteChange(self.GetReferenceValue()) * 100.0}"
        msg += f"\n\t\t violation  [%]: {abs(self.GetViolationRatio()) * 100.0}"

        return msg