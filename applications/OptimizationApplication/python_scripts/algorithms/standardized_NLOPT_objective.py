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

    def CalculateStandardizedValueAndGradients(self, control_field: numpy.ndarray, gradient_field: numpy.ndarray, save_value: bool = True) -> float:
        if self.UpdateMasterControlAndLogFields(control_field) or self.__optimization_problem.GetStep()==0:
            with TimeLogger(f"CalculateStandardizedValueAndGradients {self.GetReponse().GetName()} value", None, "Finished"):

                response_value = self.GetReponse().CalculateValue()

                if not self.__unbuffered_data.HasValue("initial_value"):
                    self.__unbuffered_data["initial_value"] = response_value

                if self.__buffered_data.HasValue("value"): del self.__buffered_data["value"]
                self.__buffered_data["value"] = response_value

                # log values
                if save_value:
                    self.LogValues()

                # print the info
                DictLogger("Objective info",self.GetInfo())

                # compute standardization factor
                self.standardization_factor = self.__scaling
                if abs(self.__unbuffered_data["initial_value"])>1e-15:
                    self.standardization_factor /= abs(self.__unbuffered_data["initial_value"])
                self.standardized_response_value = self.standardization_factor * response_value

                # compute gradients
                if gradient_field.size > 0:
                    self.gradient_collective_expression = self.CalculateGradient()
                    # save gradients
                    if save_value:
                        self.LogGradientFields()

        if gradient_field.size > 0:
            gradient_field[:] = self.standardization_factor * self.gradient_collective_expression.Evaluate().reshape(-1)

        return self.standardized_response_value

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
        if self.__optimization_problem.GetStep() > 0:
            return self.GetStandardizedValue() / self.GetInitialValue() - 1.0 if abs(self.GetInitialValue()) > 1e-12 else self.GetStandardizedValue()
        else:
            return 0.0

    def GetInfo(self) -> dict:
        info = {
            "name": self.GetReponse().GetName(),
            "type": self.__objective_type,
            "value": self.GetValue(),
            "abs_change [%]": self.GetAbsoluteChange() * 100.0,
            "rel_change [%]": self.GetRelativeChange() * 100.0
        }
        return info

    def LogValues(self) -> None:
        values_info = self.GetInfo()
        if self.__buffered_data.HasValue("rel_change [%]"): del self.__buffered_data["rel_change [%]"]
        self.__buffered_data["rel_change [%]"] = values_info["rel_change [%]"]
        if self.__buffered_data.HasValue("abs_change [%]"): del self.__buffered_data["abs_change [%]"]
        self.__buffered_data["abs_change [%]"] = values_info["abs_change [%]"]

    def LogGradientFields(self) -> None:
        # save the physical gradients for post processing in unbuffered data container.
        for physical_var, physical_gradient in self.GetRequiredPhysicalGradients().items():
            variable_name = f"d{self.GetReponse().GetName()}_d{physical_var.Name()}"
            for physical_gradient_expression in physical_gradient.GetContainerExpressions():
                if self.__unbuffered_data.HasValue(variable_name): del self.__unbuffered_data[variable_name]
                # cloning is a cheap operation, it only moves underlying pointers
                # does not create additional memory.
                self.__unbuffered_data[variable_name] = physical_gradient_expression.Clone()

        # save the filtered gradients for post processing in unbuffered data container.
        for gradient_container_expression, control in zip(self.gradient_collective_expression.GetContainerExpressions(), self.GetMasterControl().GetListOfControls()):
            variable_name = f"d{self.GetReponse().GetName()}_d{control.GetName()}"
            if self.__unbuffered_data.HasValue(variable_name): del self.__unbuffered_data[variable_name]
            # cloning is a cheap operation, it only moves underlying pointers
            # does not create additional memory.
            self.__unbuffered_data[variable_name] = gradient_container_expression.Clone()

    def LogOptimizationStep(self) -> None:
        now = datetime.datetime.now()
        now_str = now.strftime("%Y-%m-%d %H:%M:%S")
        iteration_text = f" EoF Iteration {self.__optimization_problem.GetStep()}"
        iteration_output = f"{'#'}  {iteration_text} [Time: {now_str}]  {'#'}"
        divided_line = len(iteration_output) * '#'
        to_print = f"{divided_line}\n{iteration_output}\n{divided_line}\n"
        Kratos.Logger.PrintInfo(to_print)

    def UpdateMasterControlAndLogFields(self, new_control_field: numpy.ndarray) -> bool:
        master_control_updated = False
        master_control = self.GetMasterControl()
        control_change_norm = numpy.linalg.norm(master_control.GetControlField().Evaluate().reshape(-1)-new_control_field)
        if control_change_norm > 1e-15:
            # first out put the fields and update step
            CallOnAll(self.__optimization_problem.GetListOfProcesses("output_processes"), Kratos.OutputProcess.PrintOutput)
            self.LogOptimizationStep()
            self.__optimization_problem.AdvanceStep()
            # convert numpy arry to expression
            new_control_field_exp = master_control.GetEmptyField()
            number_of_entities = []
            shapes = []
            for control in master_control.GetListOfControls():
                number_of_entities.append(len(control.GetControlField().GetContainer()))
                shapes.append(control.GetControlField().GetItemShape())
            KratosOA.CollectiveExpressionIO.Read(new_control_field_exp,new_control_field,shapes)
            # now update the master control
            master_control.Update(new_control_field_exp)
            master_control_updated = True
        return master_control_updated