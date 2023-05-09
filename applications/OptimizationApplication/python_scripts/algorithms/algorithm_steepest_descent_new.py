import KratosMultiphysics as Kratos
import KratosMultiphysics.OptimizationApplication as KratosOA
from KratosMultiphysics.python_solver import PythonSolver
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem import OptimizationProblem
from KratosMultiphysics.OptimizationApplication.algorithms.standardized_objective import StandardizedObjective
from KratosMultiphysics.OptimizationApplication.controls.master_control import MasterControl
from KratosMultiphysics.OptimizationApplication.controls.control import Control


class KratosSteepestDescent(PythonSolver):
    """
        A classical steepest descent algorithm to solve unconstrainted optimization problems.
    """

    @classmethod
    def GetDefaultParameters(cls):
        return Kratos.Parameters("""{
            "module"            : "KratosMultiphysics.OptimizationApplication.algorithms",
            "type"              : "PLEASE_PROVIDE_AN_ALGORITHM_CLASS_NAME",
            "model_part_name"   : "OptimizationModelPart",
            "objectives"        : [],
            "controls"          : [],
            "echo_level"        : 0,
            "settings"          : {
                "gradient_scaling": "inf_norm",
                "echo_level"      : 0,
                "step_size"       : float, 
                "max_iter"        : float,
            }
        }""")
    
    def __init__(self, model:Kratos.Model, parameters: Kratos.Parameters, optimization_problem: OptimizationProblem):
        super().__init__(model, parameters)

        self.model = model
        self.parameters = parameters
        self.optimization_problem = optimization_problem

        self.master_control = MasterControl() # Need to fill it with controls
        control_param_list = parameters["controls"]

        for control_param in control_param_list:
            control_name = control_param["name"].GetString()
            control = Control(control_name) # Use Control Factory 
            self.master_control.AddControl(control)

        algorithm_parameters = parameters["settings"]
        algorithm_parameters.ValidateAndAssignDefaults(self.GetDefaultParameters()["settings"])

        self.gradient_scaling = algorithm_parameters["gradient_scaling"].GetString()
        self.echo_level = algorithm_parameters["echo_level"].GetInt()

        self.step_size = algorithm_parameters["step_size"].GetInt()
        self.__max_iter = algorithm_parameters["max_iter"].GetInt()

        self.__objective = StandardizedObjective(parameters["objectives"][0], self.master_control, optimization_problem)
        self.control_field = None

    def GetMinimumBufferSize(self) -> int:
        return 2

    def Check(self):
        if len(self.parameters["objectives"]) > 1:
            raise RuntimeError(f"{self.__class__.__name__} only supports single objective optimizations.")
        
    def Initialize(self):
        super().Initialize()
        self.master_control.Initialize()

        self.converged = False
        self.control_field = self.master_control.GetEmptyControlFields() # GetInitialControlFields() later

    def Finalize(self):
        self.opt_iter = 1
        return super().Finalize()
    
    def ComputeSearchDirection(obj_grad) -> KratosOA.ContainerExpression.CollectiveExpressions:
        return obj_grad * -1.0
    
    def LineSearch(self, search_direction) -> float:
        return self.step_size / KratosOA.ContainerExpressionUtils.NormInf(search_direction)

    def SolveSolutionStep(self) -> bool:
        self.Initialize()

        while not self.converged:

            obj_val = self.__objective.CalculateStandardizedValue(self.control_field) # obj_val is typically not used. It is needed if we use Line Search Technique and for output
            obj_grad = self.__objective.CalculateStandardizedGradient()
            search_direction = self.ComputeSearchDirection(obj_grad)
            alpha = self.LineSearch(search_direction)
            self.control_field += search_direction * alpha

            self.converged = self.CheckConvergence()

        self.Finalize()

    def CheckConvergence(self) -> bool:
        return True if self.opt_iter >= self.__max_iter else False
