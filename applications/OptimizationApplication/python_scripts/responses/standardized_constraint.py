import KratosMultiphysics as Kratos
import KratosMultiphysics.OptimizationApplication as KratosOA
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem import OptimizationProblem
from KratosMultiphysics.OptimizationApplication.utilities.union_utilities import SupportedSensitivityFieldVariableTypes
from KratosMultiphysics.OptimizationApplication.responses.standardization_utilities import StandardizationUtilities

class StandardizedConstraint:
    """Standardized constraint response function

    This class creates instances to standardize any response function for the specified type of the contraint.
    Supported contraint types:
        "=",
        "<",
        "<=,
        ">",
        ">="

    The reference value for the constraint either can be the "initial_value" or a specified value.

    """
    def __init__(self, parameters: Kratos.Parameters, optimization_info: OptimizationProblem):
        default_parameters = Kratos.Parameters("""{
            "response_name": "",
            "type"         : "",
            "ref_value"    : "initial_value"
        }""")

        if parameters.Has("ref_value") and parameters["ref_value"].IsDouble():
            default_parameters["ref_value"].SetDouble(0.0)

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

        self.__utility = StandardizationUtilities(parameters["response_name"].GetString(), optimization_info)

    def GetName(self) -> str:
        return self.__utility.GetName()

    def GetType(self) -> str:
        return self.__constraint_type

    def GetReferenceValue(self):
        if self.__reference_value is None:
            self.__reference_value = self.__utility.GetScaledValue()
        return self.__reference_value

    def GetStandardizedValue(self, step_index: int = 0) -> float:
        return self.__utility.GetScaledValue(step_index, self.__standardization_value) - self.__standardization_value * self.GetReferenceValue()

    def IsActive(self, step_index: int = 0) -> bool:
        return (self.__constraint_type == "=") or self.GetStandardizedValue(step_index) >= 0.0

    def GetViolationRatio(self, step_index: int = 0) -> float:
        return (self.GetStandardizedValue(step_index) / self.GetReferenceValue() if abs(self.GetReferenceValue()) > 1e-12 else self.GetStandardizedValue(step_index)) * self.IsActive()

    def CalculateStandardizedSensitivity(self, sensitivity_variable_collective_expression_info: 'dict[SupportedSensitivityFieldVariableTypes, KratosOA.ContainerExpression.CollectiveExpressions]'):
        self.__utility.CalculateScaledSensitivity(sensitivity_variable_collective_expression_info, self.__standardization_value)

    def UpdateConstraintData(self) -> None:
        response_problem_data = self.__utility.GetBufferedData()
        response_problem_data["rel_change"] = self.__utility.GetRelativeChange()
        response_problem_data["abs_change"] = self.__utility.GetAbsoluteChange(self.GetReferenceValue())
        response_problem_data["ref_value"] = self.GetReferenceValue()
        response_problem_data["is_active"] = self.IsActive()
        response_problem_data["type"] = self.GetType()
        response_problem_data["violation"] = abs(self.GetViolationRatio())

    def GetConstraintInfo(self) -> str:
        msg = "\tConstraint info:"
        msg += f"\n\t\t name          : {self.GetName()}"
        msg += f"\n\t\t value         : {self.__utility.GetScaledValue():0.6e}"
        msg += f"\n\t\t type          : {self.GetType()}"
        msg += f"\n\t\t ref_value     : {self.GetReferenceValue():0.6e}"
        msg += f"\n\t\t is_active     : {self.IsActive()}"
        msg += f"\n\t\t rel_change [%]: {self.__utility.GetRelativeChange() * 100.0:0.6e}"
        msg += f"\n\t\t abs_change [%]: {self.__utility.GetAbsoluteChange(self.GetReferenceValue()) * 100.0:0.6e}"
        msg += f"\n\t\t violation  [%]: {abs(self.GetViolationRatio()) * 100.0:0.6e}"
        return msg