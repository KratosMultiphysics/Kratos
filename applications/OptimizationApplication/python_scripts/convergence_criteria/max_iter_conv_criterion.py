import typing
import KratosMultiphysics as Kratos
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem import OptimizationProblem
from KratosMultiphysics.OptimizationApplication.convergence_criteria.convergence_criterion import ConvergenceCriterion
from KratosMultiphysics.OptimizationApplication.utilities.logger_utilities import time_decorator

def Factory(_: Kratos.Model, parameters: Kratos.Parameters, optimization_problem: OptimizationProblem) -> ConvergenceCriterion:
    if not parameters.Has("settings"):
        raise RuntimeError(f"MaxIterConvCriterion instantiation requires a \"settings\" in parameters [ parameters = {parameters}].")
    return MaxIterConvCriterion(parameters["settings"], optimization_problem)

class MaxIterConvCriterion(ConvergenceCriterion):
    """
    MaxIterConvCriterion is a convergence criterion for optimization problems based on the maximum number of iterations.
    """
    @classmethod
    def GetDefaultParameters(cls):
        return Kratos.Parameters("""{
            "max_iter": 0
        }""")

    def __init__(self, parameters: Kratos.Parameters, optimization_problem: OptimizationProblem):
        parameters.ValidateAndAssignDefaults(self.GetDefaultParameters())

        self.__max_iter = parameters["max_iter"].GetInt()
        self.__optimization_problem = optimization_problem

        if self.__max_iter <= 0:
            raise RuntimeError(f"The number of max iterations cannot be zero or negative.")

    def Initialize(self):
        pass

    @time_decorator()
    def IsConverged(self) -> bool:
        self.__conv = (self.__optimization_problem.GetStep() >= self.__max_iter)
        return self.__conv

    def Finalize(self):
        pass

    def GetInfo(self) -> 'list[tuple[str, typing.Union[int, float, str]]]':
        info = [
                    ("type"  , "max_iter_conv_criterion"),
                    ("iter"  , f"{self.__optimization_problem.GetStep()} of {self.__max_iter}"),
                    ("status", str("converged" if self.__conv else "not converged"))
               ]
        return info