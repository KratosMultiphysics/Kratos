import KratosMultiphysics as Kratos
import KratosMultiphysics.OptimizationApplication as KratosOA
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem import OptimizationProblem
from KratosMultiphysics.OptimizationApplication.responses.response_routine import ResponseRoutine
from KratosMultiphysics.OptimizationApplication.utilities.union_utilities import SupportedSensitivityFieldVariableTypes
from KratosMultiphysics.OptimizationApplication.controls.master_control import MasterControl

class StandardizedConstraint(ResponseRoutine):
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
    def __init__(self, parameters: Kratos.Parameters, master_control: MasterControl, optimization_info: OptimizationProblem, required_buffer_size: int = 2):
        default_parameters = Kratos.Parameters("""{
            "response_name"   : "",
            "type"            : "",
            "scaling"         : 1.0,
            "scaled_ref_value": "initial_value"
        }""")

        super().__init__(master_control, parameters["response_name"].GetString(), optimization_info, required_buffer_size)

        scaling = parameters["scaling"].GetDouble()
        if scaling < 0.0:
            raise RuntimeError(f"Scaling should be always positive [ given scale = {scaling}]")

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
        if self.__constraint_type in ["<=", "="]:
            self.__scaling = scaling
        elif self.__constraint_type in [">="]:
            self.__scaling = -scaling
        else:
            raise RuntimeError(f"Provided \"type\" = {self.__constraint_type} is not supported in constraint response functions. Followings are supported options: \n\t=\n\t<=\n\t<\n\t>=\n\t>")

    def IsEqualityType(self) -> str:
        return self.__constraint_type == "="

    def GetReferenceValue(self) -> float:
        if self.__reference_value is None:
            return self.__reference_value
        else:
            raise RuntimeError(f"Response value for {self.GetReponseName()} is not calculated yet.")    

    def GetReponseValue(self, step_index):
        # return the value
        pass

    def CalculateStandardizedValue(self, control_space_updates: KratosOA.ContainerExpression.CollectiveExpressions, save_value: bool = True) -> float:
        return self.CalculateValue(control_space_updates) * self.__scaling - self.GetReferenceValue()
        
        GetScaledValue(step_index, self.__scaling) - self.__scaling * self.GetReferenceValue()
        standardized_response_value = self.CalculateValue(control_space_updates) * self.__scaling
        if not self.__initial_response_value:
            self.__initial_response_value = standardized_response_value

        # put the saving mechanism

        return standardized_response_value

    def GetViolationAbsolute(self, step_index: int = 0) -> float:
        is_violated = int(self.__constraint_type == "=") or self.GetReponseValue(step_index) >= 0.0
        return self.GetReponseValue(step_index) * is_violated

    def GetViolationRatio(self, step_index: int = 0) -> float:
        return (self.GetViolationAbsolute(step_index) / self.GetReferenceValue() if abs(self.GetReferenceValue()) > 1e-12 else self.GetViolationAbsolute(step_index))

    def CalculateStandardizedGradient(self) -> KratosOA.ContainerExpression.CollectiveExpressions:
        return self.CalculateGradient() * self.__scaling

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