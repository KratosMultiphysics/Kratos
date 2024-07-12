import KratosMultiphysics as Kratos
import KratosMultiphysics.OptimizationApplication as KratosOA
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem import OptimizationProblem
from KratosMultiphysics.OptimizationApplication.responses.response_routine import ResponseRoutine
from KratosMultiphysics.OptimizationApplication.controls.master_control import MasterControl
from KratosMultiphysics.OptimizationApplication.utilities.component_data_view import ComponentDataView
from KratosMultiphysics.OptimizationApplication.utilities.logger_utilities import DictLogger
from KratosMultiphysics.OptimizationApplication.utilities.logger_utilities import TimeLogger

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
            "response_name"       : "",
            "type"                : "",
            "scaling"             : 1.0,
            "violation_scaling"   : 1.0,
            "scaled_ref_value"    : "initial_value"
        }""")

        if parameters.Has("scaled_ref_value") and parameters["scaled_ref_value"].IsDouble():
            default_parameters["scaled_ref_value"].SetDouble(0.0)

        parameters.ValidateAndAssignDefaults(default_parameters)

        response = optimization_problem.GetResponse(parameters["response_name"].GetString())

        super().__init__(master_control, response)

        if required_buffer_size < 2:
            raise RuntimeError(f"Standardized objective requires 2 as minimum buffer size. [ response name = {self.GetResponseName()} ]")

        component_data_view = ComponentDataView(response, optimization_problem)
        component_data_view.SetDataBuffer(required_buffer_size)

        self.__optimization_problem = optimization_problem
        self.__buffered_data = component_data_view.GetBufferedData()
        self.__unbuffered_data = component_data_view.GetUnBufferedData()

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

        self.__violation_scaling = parameters["violation_scaling"].GetDouble()
        self.__constraint_type = parameters["type"].GetString()
        if self.__constraint_type in ["<=", "="]:
            self.__scaling = scaling
        elif self.__constraint_type in [">="]:
            self.__scaling = -scaling
            if not self.__reference_value is None: self.__reference_value *= -1
        else:
            raise RuntimeError(f"Provided \"type\" = {self.__constraint_type} is not supported in constraint response functions. Followings are supported options: \n\t=\n\t<=\n\t<\n\t>=\n\t>")

    def IsEqualityType(self) -> str:
        return self.__constraint_type == "="

    def GetStandardizedReferenceValue(self) -> float:
        if self.__reference_value is not None:
            return self.__reference_value
        else:
            if self.__unbuffered_data.HasValue("initial_value"):
                return self.__unbuffered_data["initial_value"] * self.__scaling
            else:
                raise RuntimeError(f"Response value for {self.GetResponseName()} is not calculated yet.")

    def CalculateStandardizedValue(self, control_field: KratosOA.CollectiveExpression, save_value: bool = True) -> float:
        with TimeLogger(f"StandardizedConstraint::Calculate {self.GetResponseName()} value", None, "Finished"):
            response_value = self.CalculateValue(control_field)
            standardized_response_value = response_value * self.__scaling

            if not self.__unbuffered_data.HasValue("initial_value"):
                self.__unbuffered_data["initial_value"] = response_value

            standardized_response_value -= self.GetStandardizedReferenceValue()

            if save_value:
                if self.__buffered_data.HasValue("value"): del self.__buffered_data["value"]
                self.__buffered_data["value"] = response_value

            DictLogger("Constraint info",self.GetInfo())

        return standardized_response_value

    def CalculateStandardizedGradient(self, save_field: bool = True) -> KratosOA.CollectiveExpression:

        with TimeLogger(f"StandardizedConstraint::Calculate {self.GetResponseName()} gradients", None, "Finished"):
            gradient_collective_expression = self.CalculateGradient()
            if save_field:
                # save the physical gradients for post processing in unbuffered data container.
                for physical_var, physical_gradient in self.GetRequiredPhysicalGradients().items():
                    variable_name = f"d{self.GetResponseName()}_d{physical_var.Name()}"
                    for physical_gradient_expression in physical_gradient.GetContainerExpressions():
                        if self.__unbuffered_data.HasValue(variable_name): del self.__unbuffered_data[variable_name]
                        # cloning is a cheap operation, it only moves underlying pointers
                        # does not create additional memory.
                        self.__unbuffered_data[variable_name] = physical_gradient_expression.Clone()

                # save the filtered gradients for post processing in unbuffered data container.
                for gradient_container_expression, control in zip(gradient_collective_expression.GetContainerExpressions(), self.GetMasterControl().GetListOfControls()):
                    variable_name = f"d{self.GetResponseName()}_d{control.GetName()}"
                    if self.__unbuffered_data.HasValue(variable_name): del self.__unbuffered_data[variable_name]
                    # cloning is a cheap operation, it only moves underlying pointers
                    # does not create additional memory.
                    self.__unbuffered_data[variable_name] = gradient_container_expression.Clone()

        return gradient_collective_expression * self.__scaling

    def GetValue(self, step_index: int = 0) -> float:
        return self.__buffered_data.GetValue("value", step_index)

    def GetStandardizedValue(self, step_index: int = 0) -> float:
        return self.GetValue(step_index) * self.__scaling - self.GetStandardizedReferenceValue()

    def GetScaledViolationValue(self, step_index: int = 0) -> float:
        value = self.GetStandardizedValue(step_index)
        return max(0.0, value * self.__violation_scaling)

    def GetAbsoluteViolation(self, step_index: int = 0) -> float:
        is_violated = self.IsEqualityType() or self.GetStandardizedValue(step_index) >= 0.0
        return self.GetStandardizedValue(step_index) * is_violated

    def GetRelativeViolation(self, step_index: int = 0) -> float:
        return self.GetAbsoluteViolation(step_index) / self.GetStandardizedReferenceValue() if abs(self.GetStandardizedReferenceValue()) > 1e-12 else self.GetAbsoluteViolation(step_index)

    def GetRelativeChange(self) -> float:
        if self.__optimization_problem.GetStep() > 1:
            return self.GetStandardizedValue() / self.GetStandardizedValue(1) - 1.0 if abs(self.GetStandardizedValue(1)) > 1e-12 else self.GetStandardizedValue()
        else:
            return 0.0

    def GetAbsoluteChange(self) -> float:
        return self.GetStandardizedValue() / self.GetStandardizedReferenceValue() - 1.0 if abs(self.GetStandardizedReferenceValue()) > 1e-12 else self.GetStandardizedValue()

    def GetInfo(self) -> dict:
        info = {
            "name": self.GetResponseName(),
            "value": self.GetValue(),
            "scaled_value": self.GetValue() * abs(self.__scaling),
            "type": self.__constraint_type,
            "ref_value": self.GetStandardizedReferenceValue(),
            "abs_change": self.GetAbsoluteChange(),
            "rel_change [%]": self.GetRelativeChange() * 100.0,
            "abs_violation [%]": self.GetAbsoluteViolation() * 100,
            "rel_violation [%]": self.GetRelativeViolation() * 100.0
        }
        return info
