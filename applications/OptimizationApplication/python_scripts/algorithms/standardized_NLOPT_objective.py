import KratosMultiphysics as Kratos
import KratosMultiphysics.OptimizationApplication as KratosOA
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem import OptimizationProblem
from KratosMultiphysics.OptimizationApplication.responses.response_routine import ResponseRoutine
from KratosMultiphysics.OptimizationApplication.controls.master_control import MasterControl
from KratosMultiphysics.OptimizationApplication.utilities.component_data_view import ComponentDataView
from KratosMultiphysics.OptimizationApplication.utilities.logger_utilities import DictLogger
from KratosMultiphysics.OptimizationApplication.utilities.logger_utilities import TimeLogger
from KratosMultiphysics.OptimizationApplication.utilities.helper_utilities import CallOnAll
import numpy
import time as timer
import datetime

class StandardizedNLOPTObjective(ResponseRoutine):
    """Standardized objective response function

    This class transformed a user-given optimization problem into the standard format.
    Supported objective types:
        "minimization",
        "maximization"

    """
    def __init__(self, parameters: Kratos.Parameters, master_control: MasterControl, optimization_problem: OptimizationProblem, required_buffer_size: int = 2):
        default_parameters = Kratos.Parameters("""{
            "response_name": "",
            "type"         : "",
            "scaling"      : 1.0
        }""")
        parameters.ValidateAndAssignDefaults(default_parameters)

        response = optimization_problem.GetResponse(parameters["response_name"].GetString())

        super().__init__(master_control, response)

        if required_buffer_size < 2:
            raise RuntimeError(f"Standardized objective requires 2 as minimum buffer size. [ response name = {self.GetReponse().GetName()} ]")

        component_data_view = ComponentDataView(response, optimization_problem)
        component_data_view.SetDataBuffer(required_buffer_size)

        self.__optimization_problem = optimization_problem
        self.__buffered_data = component_data_view.GetBufferedData()
        self.__unbuffered_data = component_data_view.GetUnBufferedData()

        scaling = parameters["scaling"].GetDouble()
        if scaling < 0.0:
            raise RuntimeError(f"Scaling should be always positive [ given scale = {scaling}]")

        self.__objective_type = parameters["type"].GetString()
        if self.__objective_type == "minimization":
            self.__scaling = scaling
        elif self.__objective_type == "maximization":
            self.__scaling = -scaling
        else:
            raise RuntimeError(f"Requesting unsupported type {self.__objective_type} for objective response function. Supported types are: \n\tminimization\n\tmaximization")

    def GetInitialValue(self) -> float:
        if self.__unbuffered_data.HasValue("initial_value"):
            return self.__unbuffered_data["initial_value"] * self.__scaling
        else:
            raise RuntimeError(f"Response value for {self.GetReponse().GetName()} is not calculated yet.")

    def Initialize(self):
        super().Initialize()

    def CalculateStandardizedValueAndGradients(self, control_field: numpy.ndarray, gradient_field: numpy.ndarray, save_value: bool = True) -> float:

        # update the master control
        update_status = self.GetMasterControl().Update(control_field)
        master_control_updated = False
        for is_updated in update_status.values():
            if is_updated:
                master_control_updated = True
                break
        # output
        if master_control_updated:
            CallOnAll(self.__optimization_problem.GetListOfProcesses("output_processes"), Kratos.OutputProcess.PrintOutput)
            self.LogOptimizationStep()
            self.__optimization_problem.AdvanceStep()

        with TimeLogger(f"CalculateStandardizedValueAndGradients {self.GetReponse().GetName()} value", None, "Finished"):
            response_value = self.CalculateValue()

            if not self.__unbuffered_data.HasValue("initial_value"):
                self.__unbuffered_data["initial_value"] = response_value

            if save_value:
                if self.__buffered_data.HasValue("value"): del self.__buffered_data["value"]
                self.__buffered_data["value"] = response_value

            # print the info
            DictLogger("Objective info",self.GetInfo())

            # compute standardization factor
            standardization_factor = self.__scaling / abs(self.__unbuffered_data["initial_value"])
            standardized_response_value = standardization_factor * response_value

            if gradient_field.size > 0:
                gradient_collective_expression = self.CalculateGradient()
                gradient_field[:] = standardization_factor * gradient_collective_expression.Evaluate()

                if save_value:
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

        return standardized_response_value

    def GetValue(self, step_index: int = 0) -> float:
        return self.__buffered_data.GetValue("value", step_index)

    def GetStandardizedValue(self, step_index: int = 0) -> float:
        return self.GetValue(step_index) * self.__scaling

    def GetRelativeChange(self) -> float:
        if self.__optimization_problem.GetStep() > 0:
            return self.GetStandardizedValue() / self.GetStandardizedValue(1) - 1.0 if abs(self.GetStandardizedValue(1)) > 1e-12 else self.GetStandardizedValue()
        else:
            return 0.0

    def GetAbsoluteChange(self) -> float:
        return self.GetValue() - self.GetInitialValue()

    def GetInfo(self) -> dict:
        info = {
            "name": self.GetReponse().GetName(),
            "type": self.__objective_type,
            "value": self.GetValue(),
            "abs_change": self.GetAbsoluteChange(),
            "rel_change [%]": self.GetRelativeChange() * 100.0
        }
        init_value = self.GetInitialValue()
        if init_value:
            info["abs_change [%]"] = self.GetAbsoluteChange()/init_value * 100

        return info

    def LogOptimizationStep(self):

        now = datetime.datetime.now()
        now_str = now.strftime("%Y-%m-%d %H:%M:%S")
        iteration_text = f" EoF Iteration {self.__optimization_problem.GetStep()}"
        iteration_output = f"{'#'}  {iteration_text} [Time: {now_str}]  {'#'}"

        divided_line = len(iteration_output) * '#'

        to_print = f"{divided_line}\n{iteration_output}\n{divided_line}\n"

        Kratos.Logger.PrintInfo(to_print)