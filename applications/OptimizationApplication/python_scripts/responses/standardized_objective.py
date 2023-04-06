import KratosMultiphysics as Kratos
import KratosMultiphysics.OptimizationApplication as KratosOA
from KratosMultiphysics.OptimizationApplication.utilities.optimization_info import OptimizationInfo
from KratosMultiphysics.OptimizationApplication.responses.response_function_data_retriever import ResponseFunctionDataRetriever
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
        self.__response_function_data_retriever = ResponseFunctionDataRetriever(self.GetName(), optimization_info)

        self.__initial_response_value = None
        self.__optimization_info = optimization_info

    def GetName(self) -> str:
        return self.__name

    def GetInitialValue(self):
        if self.__initial_response_value is None:
            self.__initial_response_value = self.GetValue()
        return self.__initial_response_value

    def GetResponseType(self) -> str:
        return self.__objective_type

    def GetStandardizedValue(self, step_index: int = 0) -> float:
        return self.__response_function_data_retriever.GetScaledValue(step_index, self.__scaling)

    def CalculateStandardizedSensitivity(self, sensitivity_variable_collective_expression_info: 'dict[SupportedSensitivityFieldVariableTypes, KratosOA.ContainerExpression.CollectiveExpressions]'):
        self.__response_function_data_retriever.CalculateScaledSensitivity(sensitivity_variable_collective_expression_info, self.__scaling)

    def UpdateOptimizationInfo(self) -> None:
        self.__optimization_info.SetValue(f"{self.__response_function_data_retriever.GetPrefix()}/type", self.GetResponseType())
        self.__optimization_info.SetValue(f"{self.__response_function_data_retriever.GetPrefix()}/relative_change", self.__response_function_data_retriever.GetRelativeChange())
        self.__optimization_info.SetValue(f"{self.__response_function_data_retriever.GetPrefix()}/absolute_change", self.__response_function_data_retriever.GetAbsoluteChange(self.GetInitialValue()))

    def GetResponseInfo(self) -> str:
        msg = "\tObjective info:"
        msg += f"\n\t\t name          : {self.GetName()}"
        msg += f"\n\t\t type          : {self.GetResponseType()}"
        msg += f"\n\t\t value         : {self.__response_function_data_retriever.GetScaledValue()}"
        msg += f"\n\t\t rel_change [%]: {self.__response_function_data_retriever.GetRelativeChange() * 100.0}"
        msg += f"\n\t\t abs_change [%]: {self.__response_function_data_retriever.GetAbsoluteChange(self.GetInitialValue()) * 100.0}"

        return msg