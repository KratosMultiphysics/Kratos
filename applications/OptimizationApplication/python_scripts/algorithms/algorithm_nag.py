import KratosMultiphysics as Kratos
import KratosMultiphysics.OptimizationApplication as KratosOA
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem import OptimizationProblem
from KratosMultiphysics.OptimizationApplication.algorithms.standardized_objective import StandardizedObjective
from KratosMultiphysics.OptimizationApplication.controls.master_control import MasterControl
from KratosMultiphysics.OptimizationApplication.algorithms.algorithm import Algorithm
from KratosMultiphysics.OptimizationApplication.utilities.component_data_view import ComponentDataView
from KratosMultiphysics.OptimizationApplication.utilities.opt_convergence import CreateConvergenceCriteria
from KratosMultiphysics.OptimizationApplication.utilities.opt_line_search import CreateLineSearch
from KratosMultiphysics.OptimizationApplication.utilities.logger_utilities import time_decorator
from KratosMultiphysics.OptimizationApplication.utilities.logger_utilities import OptimizationAlgorithmTimeLogger


def Factory(model: Kratos.Model, parameters: Kratos.Parameters, optimization_problem: OptimizationProblem):
    return AlgorithmNAG(model, parameters, optimization_problem)

class AlgorithmNAG(Algorithm):
    """
        A classical steepest descent algorithm to solve unconstrainted optimization problems.
    """

    @classmethod
    def GetDefaultParameters(cls):
        return Kratos.Parameters("""{
            "module"            : "KratosMultiphysics.OptimizationApplication.algorithms",
            "type"              : "PLEASE_PROVIDE_AN_ALGORITHM_CLASS_NAME",
            "objective"         : {},
            "controls"          : [],
            "echo_level"        : 0,
            "settings"          : {
                "echo_level"      : 0,
                "line_search"     : {},
                "conv_settings"   : {}
            }
        }""")

    def __init__(self, model:Kratos.Model, parameters: Kratos.Parameters, optimization_problem: OptimizationProblem):
        self.model = model
        self.parameters = parameters
        self._optimization_problem = optimization_problem

        parameters.ValidateAndAssignDefaults(self.GetDefaultParameters())

        self.master_control = MasterControl() # Need to fill it with controls

        for control_name in parameters["controls"].GetStringArray():
            control = optimization_problem.GetControl(control_name)
            self.master_control.AddControl(control)


        settings = parameters["settings"]
        settings.ValidateAndAssignDefaults(self.GetDefaultParameters()["settings"])

        self.echo_level = settings["echo_level"].GetInt()

        ComponentDataView("algorithm", self._optimization_problem).SetDataBuffer(self.GetMinimumBufferSize())

        self.__convergence_criteria = CreateConvergenceCriteria(settings["conv_settings"], self._optimization_problem)
        self.__line_search_method = CreateLineSearch(settings["line_search"], self._optimization_problem)

        self.__objective = StandardizedObjective(parameters["objective"], self.master_control, self._optimization_problem)
        self.__control_field = None
        self.__obj_val = None
        self.prev_update = None

        self.nu = 0.95

    def GetMinimumBufferSize(self) -> int:
        return 2

    def Check(self):
        self.master_control.Check()
        self.__objective.Check()

    @time_decorator()
    def Initialize(self):
        self.converged = False
        self.__obj_val = None
        self.master_control.Initialize()
        self.__objective.Initialize()
        self.__control_field = self.master_control.GetControlField()
        self.algorithm_data = ComponentDataView("algorithm", self._optimization_problem)

    def Finalize(self):
        self.master_control.Finalize()
        self.__objective.Finalize()

    @time_decorator()
    def ComputeSearchDirection(self, obj_grad) -> KratosOA.CollectiveExpression:
        # if self.algorithm_data.GetBufferedData().HasValue("search_direction"):
        #     search_direction = self.algorithm_data.GetBufferedData()["search_direction"] - 1.0 * obj_grad
        # else:
        #     search_direction = - obj_grad
        search_direction = - obj_grad
        self.algorithm_data.GetBufferedData()["search_direction"] = search_direction.Clone()
        return search_direction

    @time_decorator()
    def ComputeControlUpdate(self, alpha):
        search_direction = self.algorithm_data.GetBufferedData()["search_direction"]
        if isinstance(alpha, float):
            update = search_direction * alpha
        elif isinstance(alpha, KratosOA.CollectiveExpression):
            update = search_direction.Scale(alpha)
        if self.prev_update:
            mom_update = update + self.prev_update * 0.3
            self.algorithm_data.GetBufferedData()["control_field_update"] = update + mom_update * 0.3
            self.prev_update = mom_update
        else:
            self.algorithm_data.GetBufferedData()["control_field_update"] = update * 1.3
            self.prev_update = update

        for expression in self.prev_update.GetContainerExpressions():
            expression.SetExpression(expression.Flatten().GetExpression())

    @time_decorator()
    def UpdateControl(self) -> KratosOA.CollectiveExpression:
        update = self.algorithm_data.GetBufferedData()["control_field_update"]
        self.__control_field += update

    @time_decorator()
    def Output(self) -> KratosOA.CollectiveExpression:
        self.algorithm_data.GetBufferedData()["control_field"] = self.__control_field.Clone()
        self.CallOnAllProcesses(["output_processes"], Kratos.OutputProcess.PrintOutput)

    def GetCurrentObjValue(self) -> float:
        return self.__obj_val

    def GetCurrentControlField(self):
        return self.__control_field

    @time_decorator()
    def Solve(self):
        self.init_value = self.__objective.CalculateStandardizedValue(self.__control_field)
        while not self.converged:
            with OptimizationAlgorithmTimeLogger("AlgorithmNAG",self._optimization_problem.GetStep()):
                self.__obj_val = self.__objective.CalculateStandardizedValue(self.__control_field)
                obj_info = self.__objective.GetInfo()
                self.algorithm_data.GetBufferedData()["std_obj_value"] = obj_info["value"]
                self.algorithm_data.GetBufferedData()["rel_obj[%]"] = obj_info["rel_change [%]"]
                if "abs_change [%]" in obj_info:
                    self.algorithm_data.GetBufferedData()["abs_obj[%]"] = obj_info["abs_change [%]"]

                obj_grad = self.__objective.CalculateStandardizedGradient()

                self.ComputeSearchDirection(obj_grad)

                # if self.__obj_val / self.init_value < 0.01:
                #     self.__line_search_method._max_step /= 10
                #     self.init_value /= 10

                alpha = self.__line_search_method.ComputeStep()

                self.ComputeControlUpdate(alpha)

                self.Output()

                self.UpdateControl()

                self.converged = self.__convergence_criteria.IsConverged()

                self._optimization_problem.AdvanceStep()

        return self.converged

    def GetOptimizedObjectiveValue(self) -> float:
        if self.converged:
            return self.__obj_val
        else:
            raise RuntimeError("Optimization problem hasn't been solved.")
