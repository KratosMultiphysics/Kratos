import KratosMultiphysics as Kratos
import KratosMultiphysics.OptimizationApplication as KratosOA
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem import OptimizationProblem
from KratosMultiphysics.OptimizationApplication.controls.master_control import MasterControl
from KratosMultiphysics.OptimizationApplication.utilities.component_data_view import ComponentDataView
from KratosMultiphysics.OptimizationApplication.utilities.logger_utilities import DictLogger
from KratosMultiphysics.OptimizationApplication.utilities.logger_utilities import TimeLogger
from KratosMultiphysics.OptimizationApplication.utilities.response_utilities import EvaluateResponseExpression
from KratosMultiphysics.kratos_utilities import IssueDeprecationWarning
from KratosMultiphysics.OptimizationApplication.responses.response_routine import ResponseRoutine
from KratosMultiphysics.OptimizationApplication.algorithms.standardized_objective import StandardizedObjective


import numpy as np

class StandardizedSciPyObjective(ResponseRoutine):
    """Standardized SciPy objective response routine.

    This class transformed a user-given optimization problem into the standard format.
    Supported objective types:
        "minimization",
        "maximization"

    """

    def __init__(self, parameters: Kratos.Parameters, master_control: MasterControl, optimization_problem: OptimizationProblem, required_buffer_size: int = 2):
        # backward compatibility
        if parameters.Has("response_name"):
            IssueDeprecationWarning(self.__class__.__name__, "\"response_name\" is deprecated. Please use \"response_expression\".")
            parameters.AddString("response_expression", parameters["response_name"].GetString())
            parameters.RemoveValue("response_name")

        default_parameters = Kratos.Parameters("""{
            "response_expression": "",
            "type"               : "",
            "scaling"            : 1.0
        }""")
        parameters.ValidateAndAssignDefaults(default_parameters)

        response = EvaluateResponseExpression(parameters["response_expression"].GetString(), optimization_problem)

        self.objective = StandardizedObjective(parameters, master_control, optimization_problem, required_buffer_size)

        super().__init__(master_control, response)

        if required_buffer_size < 2:
            raise RuntimeError(f"Standardized objective requires 2 as minimum buffer size. [ response name = {self.GetResponseName()} ]")

        component_data_view = ComponentDataView(response, optimization_problem)
        component_data_view.SetDataBuffer(required_buffer_size)

        self.__optimization_problem = optimization_problem
        self.__buffered_data = component_data_view.GetBufferedData()
        self.__unbuffered_data = component_data_view.GetUnBufferedData()

        self.computed = False
        self.gradient_calculate_count=0

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

    def CalculateStandardizedValue(self, x:np.array, save_data: bool = True) -> float:

        control_field = self.GetMasterControl().GetEmptyField()
        shape = [c.GetItemShape() for c in control_field.GetContainerExpressions()]
        KratosOA.CollectiveExpressionIO.Read(control_field, x, shape)

        

        with TimeLogger(f"StandardizedObjective::Calculate {self.GetResponseName()} value", None, "Finished"):
            if not save_data:
                Kratos.Logger.PrintInfo("call from CalculateStandardizedGradient")
            response_value = self.objective.CalculateStandardizedValue(control_field, save_data)

            from KratosMultiphysics.OptimizationApplication.utilities.component_data_view import ComponentDataView


        
            algorithm_data = ComponentDataView("algorithm", self.__optimization_problem)
            algorithm_data.GetBufferedData()["search_direction"] = self.GetMasterControl().GetEmptyField().Clone()
            algorithm_data.GetBufferedData()["control_field_update"] = self.GetMasterControl().GetEmptyField().Clone()
            algorithm_data.GetBufferedData()["control_field"] = control_field.Clone()

            standardized_response_value = response_value * self.__scaling

            if not self.__unbuffered_data.HasValue("initial_value"):
                self.__unbuffered_data["initial_value"] = response_value
            
            if save_data:
                for process in self.__optimization_problem.GetListOfProcesses("output_processes"):
                    if process.IsOutputStep():
                        process.PrintOutput()
                self.__optimization_problem.AdvanceStep()


        return standardized_response_value

    def CalculateStandardizedGradient(self, x:np.array, save_field: bool = True) -> np.array:

        # self.CalculateStandardizedValue(x, False)  # Compute new primal if x has changed. Does nothing if x the same

        with TimeLogger(f"StandardizedObjective::Calculate {self.GetResponseName()} gradients", None, "Finished"):
            # gradient_collective_expression = self.CalculateGradient()
            gradient_collective_expression = self.objective.CalculateStandardizedGradient(save_field)
            self.gradient_calculate_count+=1
            Kratos.Logger.PrintInfo(f"Gradient expression: {self.GetResponseName()}. Gradient calcultaion count is {self.gradient_calculate_count}")
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
                    
        return gradient_collective_expression.Evaluate().reshape(-1) * self.__scaling

    def Initialize(self):
        self.objective.Initialize()

    def Check(self):
        self.objective.Check()

    def GetMasterControl(self) -> MasterControl:
        return self.objective.GetMasterControl()
                
    def Finalize(self):
        self.objective.Finalize()

    def GetStandartizedObjective(self):
        return self.objective
    
    def GetScalingFactor(self):
        return self.objective.GetScalingFactor()