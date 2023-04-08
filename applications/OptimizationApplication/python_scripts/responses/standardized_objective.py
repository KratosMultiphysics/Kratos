import KratosMultiphysics as Kratos
import KratosMultiphysics.OptimizationApplication as KratosOA
from KratosMultiphysics.OptimizationApplication.utilities.optimization_info import OptimizationInfo
from KratosMultiphysics.OptimizationApplication.utilities.union_utilities import SupportedSensitivityFieldVariableTypes
from KratosMultiphysics.OptimizationApplication.utilities.communicators.response_function_communicator import ResponseFunctionCommunicator

class StandardizedObjective:
    """Standardized objective response function

    This class creates instances to standardize any response function (make it a minimization problem).
    Supported objective types:
        "minimization",
        "maximization"

    """
    def __init__(self, parameters: Kratos.Parameters, optimization_info: OptimizationInfo):
        default_parameters = Kratos.Parameters("""{
            "response_name": "",
            "type"         : "",
            "scaling"      : 1.0
        }""")
        parameters.ValidateAndAssignDefaults(default_parameters)

        self.__objective_type = parameters["type"].GetString()
        if self.__objective_type == "minimization":
            self.__scaling = parameters["scaling"].GetDouble()
        elif self.__objective_type == "maximization":
            self.__scaling = -parameters["scaling"].GetDouble()
        else:
            raise RuntimeError(f"Requesting unsupported type {self.__objective_type} for objective response function. Supported types are: \n\tminimization\n\tmaximization")

        self.__communicator = ResponseFunctionCommunicator(parameters["response_name"].GetString(), optimization_info)
        self.__initial_response_value = None

    def GetResponseFunctionName(self) -> str:
        return self.__communicator.GetName()

    def GetInitialValue(self) -> float:
        if self.__initial_response_value is None:
            self.__initial_response_value = self.__communicator.GetScaledValue()
        return self.__initial_response_value

    def GetResponseType(self) -> str:
        return self.__objective_type

    def GetStandardizedValue(self, step_index: int = 0) -> float:
        return self.__communicator.GetScaledValue(step_index, self.__scaling)

    def CalculateStandardizedSensitivity(self, sensitivity_variable_collective_expression_info: 'dict[SupportedSensitivityFieldVariableTypes, KratosOA.ContainerExpression.CollectiveExpressions]'):
        self.__communicator.CalculateScaledSensitivity(sensitivity_variable_collective_expression_info, self.__scaling)

    def UpdateObjectiveData(self) -> None:
        response_problem_data = self.__communicator.GetBufferedDataContainer()
        response_problem_data["type"] = self.GetResponseType()
        response_problem_data["rel_change"] = self.__communicator.GetRelativeChange()
        response_problem_data["abs_change"] = self.__communicator.GetAbsoluteChange(self.GetInitialValue())

    def GetObjectiveInfo(self) -> str:
        msg = "\tObjective info:"
        msg += f"\n\t\t name          : {self.GetResponseFunctionName()}"
        msg += f"\n\t\t type          : {self.GetResponseType()}"
        msg += f"\n\t\t value         : {self.__communicator.GetScaledValue():0.6e}"
        msg += f"\n\t\t rel_change [%]: {self.__communicator.GetRelativeChange() * 100.0:0.6e}"
        msg += f"\n\t\t abs_change [%]: {self.__communicator.GetAbsoluteChange(self.GetInitialValue()) * 100.0:0.6e}"
        return msg