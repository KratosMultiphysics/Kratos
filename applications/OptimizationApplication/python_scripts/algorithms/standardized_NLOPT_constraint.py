import KratosMultiphysics as Kratos
import KratosMultiphysics.OptimizationApplication as KratosOA
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem import OptimizationProblem
from KratosMultiphysics.OptimizationApplication.responses.response_routine import ResponseRoutine
from KratosMultiphysics.OptimizationApplication.controls.master_control import MasterControl
from KratosMultiphysics.OptimizationApplication.utilities.component_data_view import ComponentDataView
from KratosMultiphysics.OptimizationApplication.utilities.logger_utilities import DictLogger
from KratosMultiphysics.OptimizationApplication.utilities.logger_utilities import TimeLogger
import numpy
import sys

class StandardizedNLOPTConstraint(ResponseRoutine):
    """Standardized constraint response function

    This class creates instances to standardize any response function for the specified type of the contraint.
    Supported contraint types:
        "=",
        "<",
        ">"

    The reference value for the constraint either can be the "initial_value" or a specified value.

    """
    def __init__(self, parameters: Kratos.Parameters, master_control: MasterControl, optimization_problem: OptimizationProblem, required_buffer_size: int = 2):
        default_parameters = Kratos.Parameters("""{
            "response_name": "",
            "type"                     : "",
            "ref_value"             : "initial_value"
        }""")

        response_name = str(parameters["response_name"].GetString())
        if not parameters.Has("ref_value"):
            raise RuntimeError(f"Constraint for response {response_name} requires reference value specified in ref_value")

        if parameters.Has("ref_value") and parameters["ref_value"].IsDouble():
            default_parameters["ref_value"].SetDouble(0.0)

        parameters.ValidateAndAssignDefaults(default_parameters)

        response = optimization_problem.GetResponse(parameters["response_name"].GetString())

        super().__init__(master_control, response)

        if required_buffer_size < 2:
            raise RuntimeError(f"Standardized constraint requires 2 as minimum buffer size. [ response name = {self.GetReponse().GetName()} ]")

        component_data_view = ComponentDataView(response, optimization_problem)
        component_data_view.SetDataBuffer(required_buffer_size)

        self.__optimization_problem = optimization_problem
        self.__buffered_data = component_data_view.GetBufferedData()
        self.__unbuffered_data = component_data_view.GetUnBufferedData()

        if parameters["ref_value"].IsDouble():
            self.__reference_value = parameters["ref_value"].GetDouble()
        elif parameters["ref_value"].IsString() and parameters["ref_value"].GetString() == "initial_value":
            self.__reference_value = None
        else:
            raise RuntimeError("Provided \"reference_type\" is not supported for constraint response functions. Followings are supported options: \n\tinitial_value\n\tfloat value")

        self.__constraint_type = parameters["type"].GetString()
        if self.__constraint_type in ["=", "<"]:
            self.__scaling = 1.0
        elif self.__constraint_type in [">"]:
            self.__scaling = -1.0
            if not self.__reference_value is None: self.__reference_value *= -1
        else:
            raise RuntimeError(f"Provided \"type\" = {self.__constraint_type} is not supported in constraint response functions. Followings are supported options: \n\t=\n\t<\n\t>")

        self.__zero_threshold = sys.float_info.epsilon

    def IsEqualityType(self) -> str:
        return self.__constraint_type == "="

    def GetReferenceValue(self) -> float:
        if self.__reference_value is not None:
            return self.__reference_value
        else:
            if self.__unbuffered_data.HasValue("initial_value"):
                self.__reference_value = self.__unbuffered_data["initial_value"]
                return self.__reference_value
            else:
                raise RuntimeError(f"Response value for {self.GetReponse().GetName()} is not calculated yet.")

    def CalculateStandardizedValueAndGradients(self, control_field: numpy.ndarray, gradient_field: numpy.ndarray, save_value: bool = True) -> float:

        with TimeLogger(f"StandardizedNLOPTConstraint::Calculate {self.GetReponse().GetName()} value", None, "Finished"):

            # first update the master control
            self.UpdateMasterControl(control_field)

            # compute the value
            response_value = self.GetReponse().CalculateValue()

            if not self.__unbuffered_data.HasValue("initial_value"):
                self.__unbuffered_data["initial_value"] = response_value

            if self.__buffered_data.HasValue("value"): del self.__buffered_data["value"]
            self.__buffered_data["value"] = response_value

            # log values
            if save_value:
                self.LogValues()

            # print the info
            DictLogger("Constraint info",self.GetInfo())

            if gradient_field.size > 0:
                gradient_collective_expression = self.CalculateGradient()
                gradient_field[:] = gradient_collective_expression.Evaluate().reshape(-1) * self.GetStandardizationFactor()

                if save_value:
                    self.LogGradientFields(gradient_collective_expression)

        return self.GetStandardizedValue()

    def GetValue(self, step_index: int = 0) -> float:
        return self.__buffered_data.GetValue("value", step_index)

    def GetStandardizationFactor(self):
        ref_value = self.GetReferenceValue()
        standardization_factor = self.__scaling
        if abs(ref_value) > self.__zero_threshold:
            standardization_factor /= abs(ref_value)
        return standardization_factor

    def GetStandardizedValue(self, step_index: int = 0) -> float:
        ref_value = self.GetReferenceValue()
        standardized_response_value = self.GetValue(step_index) - ref_value
        standardized_response_value *= self.GetStandardizationFactor()
        return standardized_response_value

    def GetAbsoluteViolation(self, step_index: int = 0) -> float:
            ref_value = self.GetReferenceValue()
            standardized_response_value = self.GetValue(step_index) - ref_value
            standardization_factor = self.__scaling
            if abs(ref_value) > self.__zero_threshold:
                standardization_factor /= abs(ref_value)
                standardization_factor *= 100
            standardized_response_value *= standardization_factor
            return standardized_response_value

    def GetInfo(self) -> dict:
        info = {
            "name": self.GetReponse().GetName(),
            "value": self.GetValue(),
            "type": self.__constraint_type,
            "ref_value": self.GetReferenceValue(),
            "violation [%]": self.GetAbsoluteViolation()
        }
        return info

    def LogValues(self) -> None:
        if self.__buffered_data.HasValue("violation [%]"): del self.__buffered_data["violation [%]"]
        self.__buffered_data["violation [%]"] = self.GetAbsoluteViolation()

    def LogGradientFields(self,gradient_collective_expression) -> None:
        # save the physical gradients for post processing in unbuffered data container.
        for physical_var, physical_gradient in self.GetRequiredPhysicalGradients().items():
            variable_name = f"d{self.GetReponse().GetName()}_d{physical_var.Name()}"
            for physical_gradient_expression in physical_gradient.GetContainerExpressions():
                if self.__unbuffered_data.HasValue(variable_name): del self.__unbuffered_data[variable_name]
                # cloning is a cheap operation, it only moves underlying pointers
                # does not create additional memory.
                self.__unbuffered_data[variable_name] = physical_gradient_expression.Clone()

        # save the filtered gradients for post processing in unbuffered data container.
        for gradient_container_expression, control in zip(gradient_collective_expression.GetContainerExpressions(), self.GetMasterControl().GetListOfControls()):
            variable_name = f"d{self.GetReponse().GetName()}_d{control.GetName()}"
            if self.__unbuffered_data.HasValue(variable_name): del self.__unbuffered_data[variable_name]
            # cloning is a cheap operation, it only moves underlying pointers
            # does not create additional memory.
            self.__unbuffered_data[variable_name] = gradient_container_expression.Clone()

    def UpdateMasterControl(self, new_control_field: numpy.ndarray) -> None:
        master_control = self.GetMasterControl()
        new_control_field_exp = master_control.GetEmptyField()
        number_of_entities = []
        shapes = []
        for control in master_control.GetListOfControls():
            number_of_entities.append(len(control.GetControlField().GetContainer()))
            shapes.append(control.GetControlField().GetItemShape())
        KratosOA.CollectiveExpressionIO.Read(new_control_field_exp,new_control_field,shapes)
        # now update the master control
        master_control.Update(new_control_field_exp)