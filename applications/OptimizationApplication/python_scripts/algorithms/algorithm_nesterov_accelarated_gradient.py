import KratosMultiphysics as Kratos
import KratosMultiphysics.OptimizationApplication as KratosOA
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem import OptimizationProblem
from KratosMultiphysics.OptimizationApplication.utilities.logger_utilities import time_decorator
from KratosMultiphysics.OptimizationApplication.algorithms.algorithm_steepest_descent import AlgorithmSteepestDescent


def Factory(model: Kratos.Model, parameters: Kratos.Parameters, optimization_problem: OptimizationProblem):
    return AlgorithmNesterovAcceleratedGradient(model, parameters, optimization_problem)

class AlgorithmNesterovAcceleratedGradient(AlgorithmSteepestDescent):
    """
        Nesterov Accelerated Gradient method to solve unconstrainted optimization problems.
        The implementation is based on https://paperswithcode.com/method/nesterov-accelerated-gradient
        In current version to save one response evaluation, the f(x) is only calculated on the momentum points.
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
                "eta"             : 0.9,
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
        # compute the correction part from momentum point
        search_direction = self.algorithm_data.GetBufferedData()["search_direction"]
        update = KratosOA.ExpressionUtils.Scale(search_direction, alpha)
        # add momentum to the correction update to compute new momentum point.
        if self.prev_update:
            mom_update = update + self.prev_update * self.eta
            self.algorithm_data.GetBufferedData()["control_field_update"] = update + mom_update * self.eta
            self.prev_update = mom_update
        else:
            self.algorithm_data.GetBufferedData()["control_field_update"] = update * (1 + self.eta)
            self.prev_update = update

        self.prev_update = KratosOA.ExpressionUtils.Collapse(self.prev_update)
