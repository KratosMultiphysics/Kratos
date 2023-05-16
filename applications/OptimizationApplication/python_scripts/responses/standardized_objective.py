import KratosMultiphysics as Kratos
import KratosMultiphysics.OptimizationApplication as KratosOA
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem import OptimizationProblem
from KratosMultiphysics.OptimizationApplication.utilities.union_utilities import SupportedSensitivityFieldVariableTypes
from KratosMultiphysics.OptimizationApplication.responses.standardization_utilities import StandardizationUtilities

class StandardizedObjective:
    """Standardized objective response function

    This class creates instances to standardize any response function (make it a minimization problem).
    Supported objective types:
        "minimization",
        "maximization"

    """
    def __init__(self, parameters: Kratos.Parameters, optimization_info: OptimizationProblem):
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

        self.__utility = StandardizationUtilities(parameters["response_name"].GetString(), optimization_info)
        self.__initial_response_value = None

    def GetName(self) -> str:
        return self.__utility.GetName()

    def GetType(self) -> str:
        return self.__objective_type

    def GetInitialValue(self) -> float:
        if self.__initial_response_value is None:
            self.__initial_response_value = self.__utility.GetScaledValue()
        return self.__initial_response_value

    def GetStandardizedValue(self, step_index: int = 0) -> float:
        """Returns the standardized objective value.

        This method returns the standardized objective value which turns all the objective types
        to minimization problems.

        If this is called with a minimization problem, then it only does the scaling.
        If this is called with a maximization problem, then it does scaling as well as negation
        to make it a minimization problem.

        Args:
            step_index (int, optional): Step index to calculate standardized value. Defaults to 0.

        Returns:
            float: Standardized resposne value.
        """
        return self.__utility.GetScaledValue(step_index, self.__scaling)

    def CalculateStandardizedSensitivity(self, sensitivity_variable_collective_expression_info: 'dict[SupportedSensitivityFieldVariableTypes, KratosOA.ContainerExpression.CollectiveExpressions]'):
        """Returns the standardized objective sensitivity.

        This method returns the standardized objective sensitivity which turns all the objective types
        to minimization problems.

        If this is called with a minimization problem, then it only does the scaling.
        If this is called with a maximization problem, then it does scaling as well as negation
        to make it a minimization problem.

        Args:
            sensitivity_variable_collective_expression_info (dict[SupportedSensitivityFieldVariableTypes, KratosOA.ContainerExpression.CollectiveExpressions]): Collective expression containing sensitivities.
        """
        self.__utility.CalculateScaledSensitivity(sensitivity_variable_collective_expression_info, self.__scaling)

    def UpdateObjectiveData(self) -> None:
        response_problem_data = self.__utility.GetBufferedData()
        response_problem_data["type"] = self.GetType()
        response_problem_data["rel_change"] = self.__utility.GetRelativeChange()
        response_problem_data["abs_change"] = self.__utility.GetAbsoluteChange(self.GetInitialValue())

    def GetObjectiveInfo(self) -> str:
        msg = "\tObjective info:"
        msg += f"\n\t\t name          : {self.GetName()}"
        msg += f"\n\t\t type          : {self.GetType()}"
        msg += f"\n\t\t value         : {self.__utility.GetScaledValue():0.6e}"
        msg += f"\n\t\t rel_change [%]: {self.__utility.GetRelativeChange() * 100.0:0.6e}"
        msg += f"\n\t\t abs_change [%]: {self.__utility.GetAbsoluteChange(self.GetInitialValue()) * 100.0:0.6e}"
        return msg