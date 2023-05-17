import numpy as np
import KratosMultiphysics as Kratos
import KratosMultiphysics.OptimizationApplication as KratosOA
from KratosMultiphysics.OptimizationApplication.controls.control import Control
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem import OptimizationProblem
from KratosMultiphysics.OptimizationApplication.algorithms.standardized_objective import StandardizedObjective
from KratosMultiphysics.OptimizationApplication.controls.master_control import MasterControl
from KratosMultiphysics.OptimizationApplication.algorithms.algorithm import Algorithm
from KratosMultiphysics.OptimizationApplication.utilities.component_data_view import ComponentDataView
from KratosMultiphysics.OptimizationApplication.utilities.helper_utilities import CallOnAll
from KratosMultiphysics import Parameters, Logger


class KratosSteepestDescent(Algorithm):
    """
        A classical steepest descent algorithm to solve unconstrainted optimization problems.
    """

    @classmethod
    def GetDefaultParameters(cls):
        return Kratos.Parameters("""{
            "module"            : "KratosMultiphysics.OptimizationApplication.algorithms",
            "type"              : "PLEASE_PROVIDE_AN_ALGORITHM_CLASS_NAME",
            "model_part_name"   : "OptimizationModelPart",
            "objective"         : {},
            "controls"          : [],
            "echo_level"        : 0,
            "settings"          : {
                "gradient_scaling": "inf_norm",
                "echo_level"      : 0,
                "step_size"       : 0.0, 
                "max_iter"        : 0
            }
        }""")
    
    def __init__(self, model:Kratos.Model, parameters: Kratos.Parameters, _optimization_problem: OptimizationProblem):
        self.model = model
        self.parameters = parameters
        self._optimization_problem = _optimization_problem

        self.master_control = MasterControl() # Need to fill it with controls
        self.__control_param_list = parameters["controls"]

        for control_param in self.__control_param_list.values():
            control_name = control_param.GetString()
            control = _optimization_problem.GetControl(control_name)
            self.master_control.AddControl(control)

        algorithm_parameters = parameters["settings"]
        algorithm_parameters.ValidateAndAssignDefaults(self.GetDefaultParameters()["settings"])

        self.gradient_scaling = algorithm_parameters["gradient_scaling"].GetString()
        self.echo_level = algorithm_parameters["echo_level"].GetInt()

        self.step_size = algorithm_parameters["step_size"].GetDouble()
        self.__max_iter = algorithm_parameters["max_iter"].GetInt()

        self.__objective = StandardizedObjective(parameters["objective"], self.master_control, self._optimization_problem)
        self.__control_field = None

        self.__obj_val = None
        self.opt_iter = None

    def GetMinimumBufferSize(self) -> int:
        return 2

    def Check(self):
        if len(self.parameters["objectives"]) > 1:
            raise RuntimeError(f"{self.__class__.__name__} only supports single objective optimizations.")
        
    def Initialize(self):
        self.converged = False
        self.__obj_val = None

        CallOnAll(self.master_control.GetListOfControls(), Control.Initialize)
        self.__objective.GetReponse().Initialize()
        self.__objective.Initialize()
        self.__objective.Check()
        self.__control_field = self.master_control.GetControlField() # GetInitialControlFields() later

    def Finalize(self):
        pass
    
    def ComputeSearchDirection(self, obj_grad) -> KratosOA.ContainerExpression.CollectiveExpressions:
        return obj_grad * -1.0
    
    def LineSearch(self, search_direction) -> float:
        norm = KratosOA.ContainerExpressionUtils.NormInf(search_direction)
        if norm:
            return self.step_size / norm
        else:
            return self.step_size
    
    def GetCurrentObjValue(self) -> float:
        return self.__obj_val
    
    def GetCurrentControlField(self):
        return self.__control_field

    def SolveOptimizationProblem(self) -> bool:
        self.Initialize()
        while not self.converged:
            print(self._optimization_problem.GetStep())

            self.__obj_val = self.__objective.CalculateStandardizedValue(self.__control_field) # __obj_val is typically not used. It is needed if we use Line Search Technique and for output
            obj_grad = self.__objective.CalculateStandardizedGradient()
            search_direction = self.ComputeSearchDirection(obj_grad)

            for control in self.master_control.GetListOfControls():
                control_data_storage = ComponentDataView(control, self._optimization_problem)
                control_data_storage.GetUnBufferedData().SetValue("search_direction", search_direction.Clone(), overwrite=True)

            alpha = self.LineSearch(search_direction)

            update = search_direction * alpha
            self.__control_field += update

            for control in self.master_control.GetListOfControls():
                control_data_storage = ComponentDataView(control, self._optimization_problem)
                control_data_storage.GetUnBufferedData().SetValue("parameter_update", update.Clone(), overwrite=True)
                control_data_storage.GetUnBufferedData().SetValue("control_field", self.__control_field.Clone(), overwrite=True)

            print(self.__objective.GetInfo())

            self.converged = self.CheckConvergence()
            self._optimization_problem.AdvanceStep()

        self.Finalize()
        return self.converged

    def CheckConvergence(self) -> bool:
        return True if self._optimization_problem.GetStep() >= self.__max_iter else False
    
    def GetOptimizedObjectiveValue(self) -> float:
        if self.converged:
            return self.__obj_val
        else:
            raise RuntimeError("Optimization problem hasn't been solved.")
