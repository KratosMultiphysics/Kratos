import KratosMultiphysics as Kratos
import KratosMultiphysics.OptimizationApplication as KratosOA
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem import OptimizationProblem
from KratosMultiphysics.OptimizationApplication.responses.response_routine import ResponseRoutine
from KratosMultiphysics.OptimizationApplication.controls.master_control import MasterControl
from KratosMultiphysics.OptimizationApplication.utilities.component_data_view import ComponentDataView
from KratosMultiphysics.OptimizationApplication.utilities.logger_utilities import DictLogger
from KratosMultiphysics.OptimizationApplication.utilities.logger_utilities import TimeLogger
from enum import Enum


class ConstraintType(Enum):
    EQUAL = "="
    LOWER = "<"
    LOWEREQ = "<="
    GREATER = ">"
    GREATEREQ = ">="

class StandardizedRGPConstraint(ResponseRoutine):
    """Standardized RGP constraint response function

    This class creates instances to standardize any response function for the specified type of the contraint.
    Supported contraint types:
        "=",
        "<",
        "<=,
        ">",
        ">="

    The reference value for the constraint either can be the "initial_value" or a specified value.

    It can be used only with relaxed_gradient_projection_algorithm. The class contains neccessary functions to compute critical zones. 

    """
    def __init__(self, parameters: Kratos.Parameters, master_control: MasterControl, optimization_problem: OptimizationProblem, required_buffer_size: int = 4):
        default_parameters = Kratos.Parameters("""{
            "response_name"       : "",
            "type"                : "",
            "scaled_ref_value"    : 0.0,
            "buffer_factor"       : 2.0,
            "max_w_c"             : 10.0,
            "tolerance"           : 1e-3
        }""")

        if parameters.Has("scaled_ref_value") and parameters["scaled_ref_value"].IsDouble():
            default_parameters["scaled_ref_value"].SetDouble(0.0)
        elif parameters.Has("scaled_ref_value") and parameters["scaled_ref_value"].GetString() == "initial_value":
            default_parameters["scaled_ref_value"].SetString("initial_value")
        else:
            type = parameters["scaled_ref_value"].GetString()
            raise RuntimeError(f"Provided \"scaled_ref_value\" = {type} is not supported for constraint response functions. Followings are supported options: \n\tinitial_value\n\tfloat value")

        parameters.ValidateAndAssignDefaults(default_parameters)

        response = optimization_problem.GetResponse(parameters["response_name"].GetString())

        super().__init__(master_control, response)

        if required_buffer_size < 4:
            raise RuntimeError(f"Standardized RGP constraint requires 4 as minimum buffer size. [ response name = {self.GetResponseName()} ]")

        component_data_view = ComponentDataView(response, optimization_problem)
        component_data_view.SetDataBuffer(required_buffer_size)

        self.__optimization_problem = optimization_problem
        self.__buffered_data = component_data_view.GetBufferedData()
        self.__unbuffered_data = component_data_view.GetUnBufferedData()

        if parameters["scaled_ref_value"].IsDouble():
            self.__reference_value = parameters["scaled_ref_value"].GetDouble()
        elif parameters["scaled_ref_value"].IsString() and parameters["scaled_ref_value"].GetString() == "initial_value":
            self.__reference_value = None


        self.__constraint_type = ConstraintType(parameters["type"].GetString())
        if self.__constraint_type in [ConstraintType.EQUAL, ConstraintType.LOWEREQ, ConstraintType.LOWER]:
            self.__scaling = 1.0
        elif self.__constraint_type in [ConstraintType.GREATEREQ, ConstraintType.GREATER]:
            self.__scaling = -1.0
            if not self.__reference_value is None: self.__reference_value *= -1
        else:
            raise RuntimeError(f"Provided \"type\" = {self.__constraint_type} is not supported in constraint response functions. Followings are supported options: \n\t=\n\t<=\n\t<\n\t>=\n\t>")

        # buffer coefficients
        self.BSF = parameters["buffer_factor"].GetDouble()
        self.BSF_init = self.BSF
        self.CBV = 0.0
        self.BS = 1e-6
        self.max_w_c = parameters["max_w_c"].GetDouble()
        self.CF = 1.0

        self.tolerance = parameters["tolerance"].GetDouble()

    def IsEqualityType(self) -> str:
        return self.__constraint_type == ConstraintType.EQUAL

    def ComputeW(self) -> float:
        if self.IsEqualityType():
            return max(abs(self.GetStandardizedValue() - (self.CBV - self.BS)), 0.0) / self.BS
        else:
            return max(self.GetStandardizedValue() - (self.CBV - self.BS), 0.0) / self.BS

    def PredictW(self, value):
        if self.IsEqualityType():
            return max(abs(value - (self.CBV - self.BS)), 0.0) / self.BS
        else:
            return max(value - (self.CBV - self.BS), 0.0) / self.BS

    def Compute_W_relax(self, w:float) -> float:
        if w <= 0.0:
            return 0.0
        elif w <= 1.0:
            return w
        else:
            return 1.0

    def Compute_W_correction(self, w:float) -> float:
        if w <= 1.0:
            return 0.0
        elif w <= self.max_w_c:
            return (w - 1.0) * self.BSF_init * self.CF
        else:
            return (self.max_w_c - 1) * self.BSF_init * self.CF

    def UpdateBufferSize(self):
        step = self.__optimization_problem.GetStep() + 1
        if step == 1:
            return

        max_step = min(step, self.__buffered_data.GetBufferSize())
        values = []
        delta_values = []
        abs_delta_values = []
        for index in range(max_step):
            values.append(self.GetStandardizedValue(index))
            if index > 0:
                delta_values.append(values[index - 1] - values[index])
                abs_delta_values.append(abs(values[index - 1] - values[index]))

        max_delta = max(abs_delta_values)

        # # Update CBV if violating
        if self.IsEqualityType():
            if step > 1:
                if values[0] > 0.0 and values[1] > 0.0 and delta_values[0] > 0.0:
                    self.CF += 0.5

                if values[0] < 0.0 and values[1] < 0.0 and delta_values[0] < 0.0:
                    self.CF += 0.5
        else:
            if step > 1:
                if values[0] > 0.0 and values[1] > 0.0 and delta_values[0] > 0.0:
                    self.CF += 0.5

                if values[0] < 0.0 and values[1] < 0.0 and delta_values[0] < 0.0:
                    self.CF -= 0.5
                    self.CF = max (self.CF, 1.0)

        self.BS = self.BSF * max_delta
        # Update BSF if zig zaging
        if step > 3:
            if delta_values[0] * delta_values[1] < 0.0 and delta_values[1] * delta_values[2] < 0.0:
                self.BSF += 1

        print(f"RGP Constraint {self.GetResponseName()}:: UpdateBufferSize")
        print(f"RGP Constraint {self.GetResponseName()}:: CBV = {self.CBV}")
        print(f"RGP Constraint {self.GetResponseName()}:: BSF = {self.BSF}")
        print(f"RGP Constraint {self.GetResponseName()}:: BS = {self.BS}")
        print(f"RGP Constraint {self.GetResponseName()}:: CF = {self.CF}")

    def IsActiveConstrant(self):
        return self.ComputeW() > 0.0

    def IsSatisfied(self, value=None):
        if not value:
            value = abs(self.GetStandardizedValue()) if self.IsEqualityType() else self.GetStandardizedValue()
        if self.IsEqualityType():
            if abs(value) <= self.tolerance: return True
            else: return False
        else:
            if value <= self.tolerance: return True
            else: return False

    def GetReferenceValue(self) -> float:
        if self.__reference_value is not None:
            return self.__reference_value
        else:
            if self.__unbuffered_data.HasValue("initial_value"):
                return self.__unbuffered_data["initial_value"] * self.__scaling
            else:
                raise RuntimeError(f"Response value for {self.GetResponseName()} is not calculated yet.")

    def CalculateStandardizedValue(self, control_field: KratosOA.CollectiveExpression, save_value: bool = True) -> float:
        with TimeLogger(f"StandardizedRGPConstraint::Calculate {self.GetResponseName()} value", None, "Finished"):
            response_value = self.CalculateValue(control_field)
            standardized_response_value = response_value * self.__scaling

            if not self.__unbuffered_data.HasValue("initial_value"):
                self.__unbuffered_data["initial_value"] = response_value

            standardized_response_value -= self.GetReferenceValue()

            if save_value:
                if self.__buffered_data.HasValue("value"): del self.__buffered_data["value"]
                self.__buffered_data["value"] = response_value

            DictLogger("Constraint info",self.GetInfo())

        return standardized_response_value

    def CalculateStandardizedGradient(self, save_field: bool = True) -> KratosOA.CollectiveExpression:

        with TimeLogger(f"StandardizedRGPConstraint::Calculate {self.GetResponseName()} gradients", None, "Finished"):
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

        if self.IsEqualityType():
            factor = 1.0 if self.GetStandardizedValue() > 0.0 else -1.0
        else: factor = 1.0
        return gradient_collective_expression * self.__scaling * factor

    def GetValue(self, step_index: int = 0) -> float:
        return self.__buffered_data.GetValue("value", step_index)

    def GetStandardizedValue(self, step_index: int = 0) -> float:
        return self.GetValue(step_index) * self.__scaling - self.GetReferenceValue()

    def GetIfConstraintActive(self, step_index: int = 0) -> bool:
        return True

    def GetAbsoluteViolation(self, step_index: int = 0) -> float:
        is_violated = self.IsEqualityType() or self.GetStandardizedValue(step_index) >= 0.0
        return self.GetStandardizedValue(step_index) * is_violated

    def GetRelativeViolation(self, step_index: int = 0) -> float:
        return self.GetAbsoluteViolation(step_index) / self.GetReferenceValue() if abs(self.GetReferenceValue()) > 1e-12 else None

    def GetRelativeChange(self) -> float:
        if self.__optimization_problem.GetStep() > 1:
            return self.GetStandardizedValue() / self.GetStandardizedValue(1) - 1.0 if abs(self.GetStandardizedValue(1)) > 1e-12 else self.GetStandardizedValue()
        else:
            return 0.0

    def GetAbsoluteChange(self) -> float:
        return self.GetValue() - self.__unbuffered_data["initial_value"] if abs(self.GetReferenceValue()) > 1e-12 else self.GetStandardizedValue()

    def GetInfo(self) -> dict:
        info = {
            "name": self.GetResponseName(),
            "value": self.GetValue(),
            "type": self.__constraint_type,
            "ref_value": self.GetReferenceValue(),
            "standartiezed_value": self.GetStandardizedValue(),
            "abs_change": self.GetAbsoluteChange(),
            "abs_violation": self.GetAbsoluteViolation(),
            "rel_change [%]": self.GetRelativeChange() * 100.0,
            "rel_violation [%]": self.GetRelativeViolation() * 100.0
        }
        return info
