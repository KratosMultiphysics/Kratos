import KratosMultiphysics as Kratos
import KratosMultiphysics.OptimizationApplication as KratosOA
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem import OptimizationProblem
from KratosMultiphysics.OptimizationApplication.responses.response_routine import ResponseRoutine
from KratosMultiphysics.OptimizationApplication.controls.master_control import MasterControl
from KratosMultiphysics.OptimizationApplication.utilities.component_data_view import ComponentDataView
from KratosMultiphysics.OptimizationApplication.utilities.logger_utilities import DictLogger
from KratosMultiphysics.OptimizationApplication.utilities.logger_utilities import TimeLogger
from KratosMultiphysics.OptimizationApplication.utilities.response_utilities import EvaluateResponseExpression
from KratosMultiphysics.kratos_utilities import IssueDeprecationWarning
import numpy as np

class StandardizedSciPyConstraint(ResponseRoutine):
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
        # backward compatibility
        if parameters.Has("response_name"):
            IssueDeprecationWarning(self.__class__.__name__, "\"response_name\" is deprecated. Please use \"response_expression\".")
            parameters.AddString("response_expression", parameters["response_name"].GetString())
            parameters.RemoveValue("response_name")

        default_parameters = Kratos.Parameters("""{
            "response_expression" : "",
            "upper_boundary"      : 1e16,
            "lower_boundary"      :-1e16
        }""")

        parameters.ValidateAndAssignDefaults(default_parameters)

        response = EvaluateResponseExpression(parameters["response_expression"].GetString(), optimization_problem)

        super().__init__(master_control, response)

        if required_buffer_size < 2:
            raise RuntimeError(f"Standardized objective requires 2 as minimum buffer size. [ response name = {self.GetResponseName()} ]")

        component_data_view = ComponentDataView(response, optimization_problem)
        component_data_view.SetDataBuffer(required_buffer_size)

        self.__optimization_problem = optimization_problem
        self.__buffered_data = component_data_view.GetBufferedData()
        self.__unbuffered_data = component_data_view.GetUnBufferedData()

        self.upper_bound = parameters["upper_boundary"].GetDouble()
        self.lower_bound = parameters["lower_boundary"].GetDouble()

    def CalculateStandardizedValue(self, x:np.array, save_value: bool = True) -> float:

        control_field = self.GetMasterControl().GetEmptyField()
        shape = [c.GetItemShape() for c in control_field.GetContainerExpressions()]
        KratosOA.CollectiveExpressionIO.Read(control_field, x, shape)

        with TimeLogger(f"StandardizedObjective::Calculate {self.GetResponseName()} value", None, "Finished"):
            response_value = self.CalculateValue(control_field)
            standardized_response_value = response_value

            if not self.__unbuffered_data.HasValue("initial_value"):
                self.__unbuffered_data["initial_value"] = response_value

            if save_value:
                if self.__buffered_data.HasValue("value"): del self.__buffered_data["value"]
                self.__buffered_data["value"] = response_value

            DictLogger("Objective info",self.GetInfo())

        return standardized_response_value

    def GetValue(self, step_index: int = 0) -> float:
        return self.__buffered_data.GetValue("value", step_index)

    def GetStandardizedValue(self, step_index: int = 0) -> float:
        return self.GetValue(step_index)

    def CalculateStandardizedGradient(self, x:np.array, save_field: bool = True) -> np.array:

        self.CalculateStandardizedValue(x, False)  # Compute new primal if x has changed. Does nothing if x the same

        with TimeLogger(f"StandardizedObjective::Calculate {self.GetResponseName()} gradients", None, "Finished"):
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
                    
        return gradient_collective_expression.Evaluate()

    def GetValue(self, step_index: int = 0) -> float:
        return self.__buffered_data.GetValue("value", step_index)

    def GetStandardizedValue(self, step_index: int = 0) -> float:
        return self.GetValue(step_index) - self.GetStandardizedReferenceValue()
    
    def GetUpperBound(self):
        return self.upper_bound
    
    def GetLowerBound(self):
        return self.lower_bound

    def GetInfo(self) -> dict:
        info = {
            "name": self.GetResponseName(),
            "value": self.GetValue(),
            "upper_bound": self.upper_bound,
            "lower_bound": self.lower_bound
        }
        return info
