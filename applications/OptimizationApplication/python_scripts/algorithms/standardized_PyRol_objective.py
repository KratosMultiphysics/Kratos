import KratosMultiphysics as Kratos
import KratosMultiphysics.OptimizationApplication as KratosOA
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem import OptimizationProblem
from KratosMultiphysics.OptimizationApplication.algorithms.standardized_objective import StandardizedObjective
from KratosMultiphysics.OptimizationApplication.controls.master_control import MasterControl
from KratosMultiphysics.OptimizationApplication.utilities.component_data_view import ComponentDataView
from KratosMultiphysics.OptimizationApplication.utilities.logger_utilities import DictLogger
from KratosMultiphysics.OptimizationApplication.utilities.logger_utilities import TimeLogger
from pyrol import Objective
import numpy as np

class StandardizedPyRolObjective(Objective):
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

        self.__objective = StandardizedObjective(parameters, master_control, optimization_problem, required_buffer_size)

        super().__init__()

    def value(self, x: np.array, tol: float, save_data: bool = True) -> float:
        control_field = self.GetMasterControl().GetEmptyField()
        shape = [c.GetItemShape() for c in control_field.GetContainerExpressions()]
        KratosOA.CollectiveExpressionIO.Read(control_field, x, shape)
        value = self.__objective.CalculateStandardizedValue(control_field, save_data)

        if save_data:
            for process in self.__optimization_problem.GetListOfProcesses("output_processes"):
                if process.IsOutputStep():
                    process.PrintOutput()

            self.__optimization_problem.AdvanceStep()
        return value

    def gradient(self, g:np.array, x:np.array, save_field: bool = True):

        self.value(x, False)  # Compute new primal if x has changed. Does nothing if x the same

        result = self.__objective.CalculateStandardizedGradient(save_field)

        g = result.Evaluate().reshape(-1) * self.__objective.GetScalingFactor()

    #     with TimeLogger(f"StandardizedObjective::Calculate {self.GetResponseName()} gradients", None, "Finished"):
    #         gradient_collective_expression = self.CalculateGradient()
    #         if save_field:
    #             # save the physical gradients for post processing in unbuffered data container.
    #             for physical_var, physical_gradient in self.GetRequiredPhysicalGradients().items():
    #                 variable_name = f"d{self.GetResponseName()}_d{physical_var.Name()}"
    #                 for physical_gradient_expression in physical_gradient.GetContainerExpressions():
    #                     if self.__unbuffered_data.HasValue(variable_name): del self.__unbuffered_data[variable_name]
    #                     # cloning is a cheap operation, it only moves underlying pointers
    #                     # does not create additional memory.
    #                     self.__unbuffered_data[variable_name] = physical_gradient_expression.Clone()

    #             # save the filtered gradients for post processing in unbuffered data container.
    #             for gradient_container_expression, control in zip(gradient_collective_expression.GetContainerExpressions(), self.GetMasterControl().GetListOfControls()):
    #                 variable_name = f"d{self.GetResponseName()}_d{control.GetName()}"
    #                 if self.__unbuffered_data.HasValue(variable_name): del self.__unbuffered_data[variable_name]
    #                 # cloning is a cheap operation, it only moves underlying pointers
    #                 # does not create additional memory.
    #                 self.__unbuffered_data[variable_name] = gradient_container_expression.Clone()

    #     g = gradient_collective_expression * self.__scaling

    # def hessVec(self, hv, v, x, tol):
    #     raise RuntimeError("Hessian-vector product is not implemented for the standardized objective response function.")

    # def GetRelativeChange(self) -> float:
    #     if self.__optimization_problem.GetStep() > 0:
    #         return self.GetStandardizedValue() / self.GetStandardizedValue(1) - 1.0 if abs(self.GetStandardizedValue(1)) > 1e-12 else self.GetStandardizedValue()
    #     else:
    #         return 0.0

    # def GetAbsoluteChange(self) -> float:
    #     return self.GetValue() - self.GetInitialValue()

    # def GetInfo(self) -> dict:
    #     info = {
    #         "name": self.GetResponseName(),
    #         "type": self.__objective_type,
    #         "value": self.GetValue(),
    #         "std_value": self.GetStandardizedValue(),
    #         "abs_change": self.GetAbsoluteChange(),
    #         "rel_change [%]": self.GetRelativeChange() * 100.0
    #     }
    #     init_value = self.GetInitialValue()
    #     if init_value:
    #         info["abs_change [%]"] = self.GetAbsoluteChange()/init_value * 100

    #     return info