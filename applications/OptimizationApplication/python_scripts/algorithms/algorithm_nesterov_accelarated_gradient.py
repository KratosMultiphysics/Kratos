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
from KratosMultiphysics.OptimizationApplication.algorithms.algorithm_steepest_descent import AlgorithmSteepestDescent


def Factory(model: Kratos.Model, parameters: Kratos.Parameters, optimization_problem: OptimizationProblem):
    return AlgorithmNesterovAcceleratedGradient(model, parameters, optimization_problem)

class AlgorithmNesterovAcceleratedGradient(AlgorithmSteepestDescent):
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
                "eta"             : 0.95,
                "echo_level"      : 0,
                "line_search"     : {},
                "conv_settings"   : {}
            }
        }""")

    def __init__(self, model:Kratos.Model, parameters: Kratos.Parameters, optimization_problem: OptimizationProblem):
        super().__init__(model, parameters, optimization_problem)
        self.prev_update = None
        self.eta = self.parameters["settings"]["eta"].GetDouble()

    @time_decorator()
    def ComputeControlUpdate(self, alpha):
        search_direction = self.algorithm_data.GetBufferedData()["search_direction"]
        if isinstance(alpha, float):
            update = search_direction * alpha
        elif isinstance(alpha, KratosOA.CollectiveExpression):
            update = search_direction.Scale(alpha)
        if self.prev_update:
            mom_update = update + self.prev_update * self.eta
            self.algorithm_data.GetBufferedData()["control_field_update"] = update + mom_update * self.eta
            self.prev_update = mom_update
        else:
            self.algorithm_data.GetBufferedData()["control_field_update"] = update * (1 + self.eta)
            self.prev_update = update

        for expression in self.prev_update.GetContainerExpressions():
            expression.SetExpression(expression.Flatten().GetExpression())
