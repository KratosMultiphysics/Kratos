import KratosMultiphysics as Kratos
import KratosMultiphysics.OptimizationApplication as KratosOA
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem import OptimizationProblem
from KratosMultiphysics.OptimizationApplication.responses.response_routine import ResponseRoutine
from KratosMultiphysics.OptimizationApplication.controls.master_control import MasterControl
from KratosMultiphysics.OptimizationApplication.utilities.component_data_view import ComponentDataView

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
    def __init__(self, parameters: Kratos.Parameters, master_control: MasterControl, optimization_problem: OptimizationProblem, required_buffer_size: int = 2):
        default_parameters = Kratos.Parameters("""{
            "response_name"   : "",
            "type"            : "",
            "scaling"         : 1.0,
            "scaled_ref_value": "initial_value"
        }""")

        if parameters.Has("ref_value") and parameters["ref_value"].IsDouble():
            default_parameters["ref_value"].SetDouble(0.0)

        parameters.ValidateAndAssignDefaults(default_parameters)

        response = optimization_problem.GetResponse(parameters["response_name"].GetString())

        super().__init__(master_control, response)

        if required_buffer_size < 2:
            raise RuntimeError(f"Standardized objective requires 2 as minimum buffer size. [ response name = {self.GetReponse().GetName()} ]")

        self.__optimization_problem = optimization_problem
        self.__component_data_view = ComponentDataView(response, optimization_problem)
        self.__component_data_view.SetDataBuffer(required_buffer_size)

        scaling = parameters["scaling"].GetDouble()
        if scaling < 0.0:
            raise RuntimeError(f"Scaling should be always positive [ given scale = {scaling}]")

        if parameters["scaled_ref_value"].IsDouble():
            self.__ref_type = "specified_value"
            self.__reference_value = parameters["scaled_ref_value"].GetDouble()
        elif parameters["scaled_ref_value"].IsString() and parameters["scaled_ref_value"].GetString() == "initial_value":
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
        if self.__reference_value is not None:
            return self.__reference_value
        else:
            raise RuntimeError(f"Response value for {self.GetReponse().GetName()} is not calculated yet.")

    def CalculateStandardizedValue(self, control_field: KratosOA.ContainerExpression.CollectiveExpressions, save_value: bool = True) -> float:
        response_value = self.CalculateValue(control_field)
        standardized_response_value = response_value * self.__scaling

        if self.__reference_value is None:
            self.__reference_value = standardized_response_value

        standardized_response_value -= self.GetReferenceValue()

        if save_value:
            data = self.__component_data_view.GetBufferedData()
            data["value"] = response_value
            data["standardized_value"] = standardized_response_value

        return standardized_response_value

    def CalculateStandardizedGradient(self, save_value: bool = True) -> KratosOA.ContainerExpression.CollectiveExpressions:
        gradient_collective_expression = self.CalculateGradient()

        if save_value:
            for gradient_container_expression, control in zip(gradient_collective_expression.GetContainerExpressions(), self.GetMasterControl().GetListOfControls()):
                self.__component_data_view.GetUnBufferedData()[f"d{self.GetReponse().GetName()}_d{control.GetName()}"] = gradient_container_expression.Clone()

        return gradient_collective_expression * self.__scaling

    def GetValue(self, step_index: int) -> float:
        return self.__component_data_view.GetBufferedData().GetValue("value", step_index)

    def GetStandardizedValue(self, step_index: int) -> float:
        return self.__component_data_view.GetBufferedData().GetValue("standardized_value", step_index)

    def GetAbsoluteViolation(self, step_index: int = 0) -> float:
        is_violated = self.IsEqualityType() or self.GetStandardizedValue(step_index) >= 0.0
        return self.GetStandardizedValue(step_index) * is_violated

    def GetRelativeViolation(self, step_index: int = 0) -> float:
        is_violated = self.IsEqualityType() or self.GetStandardizedValue(step_index) >= 0.0
        return ((self.GetAbsoluteViolation(step_index) / self.GetReferenceValue() if abs(self.GetReferenceValue()) > 1e-12 else self.GetAbsoluteViolation(step_index))) * is_violated

    def GetRelativeChange(self) -> float:
        if self.__optimization_problem.GetStep() > 1:
            return self.GetStandardizedValue() / self.GetStandardizedValue(1) - 1.0 if abs(self.GetStandardizedValue(1)) > 1e-12 else self.GetStandardizedValue()
        else:
            return 0.0

    def GetAbsoluteChange(self) -> float:
        return self.GetStandardizedValue() / self.GetReferenceValue() - 1.0 if abs(self.GetReferenceValue()) > 1e-12 else self.GetStandardizedValue()

    def UpdateOptimizationProblemData(self) -> None:
        response_problem_data = self.__component_data_view.GetBufferedData()
        response_problem_data["type"] = self.__constraint_type
        response_problem_data["ref_value"] = self.GetReferenceValue()
        response_problem_data["rel_change"] = self.GetRelativeChange()
        response_problem_data["abs_change"] = self.GetAbsoluteChange()
        response_problem_data["rel_violation"] = self.GetRelativeViolation()
        response_problem_data["abs_violation"] = self.GetAbsoluteViolation()

    def GetInfo(self) -> str:
        msg = "\tConstraint info:"
        msg += f"\n\t\t name             : {self.GetReponse().GetName()}"
        msg += f"\n\t\t value            : {self.GetValue():0.6e}"
        msg += f"\n\t\t type             : {self.__constraint_type}"
        msg += f"\n\t\t ref_value        : {self.GetReferenceValue():0.6e}"
        msg += f"\n\t\t abs_change       : {self.GetAbsoluteChange():0.6e}"
        msg += f"\n\t\t rel_change [%]   : {self.GetRelativeChange() * 100.0:0.6e}"
        msg += f"\n\t\t abs_violation    : {self.GetAbsoluteViolation():0.6e}"
        msg += f"\n\t\t rel_violation [%]: {self.GetRelativeViolation() * 100.0:0.6e}"
        return msg