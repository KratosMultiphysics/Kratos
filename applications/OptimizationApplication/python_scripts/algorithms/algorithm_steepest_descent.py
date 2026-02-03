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
    return AlgorithmSteepestDescent(model, parameters, optimization_problem)

class AlgorithmSteepestDescent(Algorithm):
    """
        A classical steepest descent algorithm to solve unconstrained optimization problems.
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
        self._optimization_problem.AddComponent(self.master_control)

        for control_name in parameters["controls"].GetStringArray():
            control = optimization_problem.GetControl(control_name)
            ComponentDataView(control, self._optimization_problem).SetDataBuffer(self.GetMinimumBufferSize())
            self.master_control.AddControl(control)


        settings = parameters["settings"]
        settings.ValidateAndAssignDefaults(self.GetDefaultParameters()["settings"])

        self.echo_level = settings["echo_level"].GetInt()

        ComponentDataView("algorithm", self._optimization_problem).SetDataBuffer(self.GetMinimumBufferSize())

        self.__convergence_criteria = CreateConvergenceCriteria(settings["conv_settings"], self._optimization_problem)
        self.__line_search_method = CreateLineSearch(settings["line_search"], self._optimization_problem)

        self.__objective = StandardizedObjective(parameters["objective"], self.master_control, self._optimization_problem)
        self._optimization_problem.AddComponent(self.__objective)
        self.__control_field = None
        self.__obj_val = None

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
        search_direction = obj_grad * -1.0

        E_dependant = True
        if E_dependant:

            print("Before")
            print(search_direction.GetContainerExpressions()[0].Evaluate())
            print(search_direction.GetContainerExpressions()[1].Evaluate())

            E_x_container_exp = search_direction.GetContainerExpressions()[0].Clone()
            E_y_container_exp = search_direction.GetContainerExpressions()[1].Clone()

            print(E_x_container_exp)
            print(E_y_container_exp)

            avg_exp = Kratos.Expression.Utils.Scale((E_x_container_exp + E_y_container_exp), 0.5).Clone()
            
            E_x_container_exp_orig = search_direction.GetContainerExpressions()[0]
            E_y_container_exp_orig = search_direction.GetContainerExpressions()[1]

            E_x_container_exp_orig.SetExpression(Kratos.Expression.Utils.Collapse(avg_exp).GetExpression())
            E_y_container_exp_orig.SetExpression(Kratos.Expression.Utils.Collapse(avg_exp).GetExpression())
            # for container_expression in search_direction.GetContainerExpressions():
            #     current_gradient : Kratos.Expression = container_expression.Clone()
            #     print(current_gradient)
            #     current_gradient.
            print("After")
            print(search_direction.GetContainerExpressions()[0].Evaluate())
            print(search_direction.GetContainerExpressions()[1].Evaluate())


        #raise RuntimeError(222222)

        self.master_control.EmplaceData("search_direction", search_direction, self._optimization_problem)
        return search_direction

    @time_decorator()
    def ComputeControlUpdate(self, alpha):
        search_direction = self.master_control.GetData("search_direction", self._optimization_problem)
        update = KratosOA.ExpressionUtils.Scale(search_direction, alpha)
        print("searchdir_new:\n", search_direction)
        print("update_vor_search:\n", update)
        print("alpha:\n", alpha)
        self.master_control.EmplaceData("control_field_update", update, self._optimization_problem)

    @time_decorator()
    def UpdateControl(self) -> KratosOA.CollectiveExpression:
        update = self.master_control.GetData("control_field_update", self._optimization_problem)
        print("control_field:\n", self.__control_field)
        print("update:\n", update)
        self.__control_field = KratosOA.ExpressionUtils.Collapse(self.__control_field + update)

    @time_decorator()
    def Output(self) -> KratosOA.CollectiveExpression:
        self.master_control.EmplaceData("control_field", self.__control_field.Clone(), self._optimization_problem)
        for process in self._optimization_problem.GetListOfProcesses("output_processes"):
            if process.IsOutputStep():
                process.PrintOutput()

    def GetCurrentObjValue(self) -> float:
        return self.__obj_val

    def GetCurrentControlField(self):
        return self.__control_field

    @time_decorator()
    def Solve(self):
        while not self.converged:
            with OptimizationAlgorithmTimeLogger("AlgorithmSteepestDescent",self._optimization_problem.GetStep()):
                self._InitializeIteration()

                self.__obj_val = self.__objective.CalculateStandardizedValue(self.__control_field)
                obj_info = self.__objective.GetInfo()
                self.algorithm_data.GetBufferedData()["std_obj_value"] = obj_info["std_value"]
                self.algorithm_data.GetBufferedData()["rel_obj[%]"] = obj_info["rel_change [%]"]
                if "abs_change [%]" in obj_info:
                    self.algorithm_data.GetBufferedData()["abs_obj[%]"] = obj_info["abs_change [%]"]

                obj_grad = self.__objective.CalculateStandardizedGradient()

                self.ComputeSearchDirection(obj_grad)

                alpha = self.__line_search_method.ComputeStep(self.master_control)

                self.ComputeControlUpdate(alpha)

                self._FinalizeIteration()

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
