import KratosMultiphysics as Kratos
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem import OptimizationProblem
from KratosMultiphysics.OptimizationApplication.utilities.component_data_view import ComponentDataView
import KratosMultiphysics.OptimizationApplication as KratosOA
from KratosMultiphysics.OptimizationApplication.utilities.logger_utilities import DictLogger
from KratosMultiphysics.OptimizationApplication.utilities.logger_utilities import TimeLogger

def CreateConvergenceCriteria(parameters: Kratos.Parameters, optimization_problem: OptimizationProblem):
    type = parameters["type"].GetString()
    if type == "max_iter":
        return MaxIterConvCriterion(parameters, optimization_problem)
    elif type == "l2_norm":
        return L2ConvCriterion(parameters, optimization_problem)
    else:
        raise RuntimeError(f"CreateConvergenceCriteria: unsupported convergence type {type}.")

class MaxIterConvCriterion:
    @classmethod
    def GetDefaultParameters(cls):
        return Kratos.Parameters("""{
            "type"              : "max_iter",
            "max_iter"          : 0
        }""")

    def __init__(self, parameters: Kratos.Parameters, optimization_problem: OptimizationProblem):
        parameters.ValidateAndAssignDefaults(self.GetDefaultParameters())
        self.__max_iter = parameters["max_iter"].GetInt()
        self.__optimization_problem = optimization_problem

    def IsConverged(self, search_direction=None) -> bool:
        with TimeLogger("MaxIterConvCriterion::IsConverged", None, "Finished"):
            iter = self.__optimization_problem.GetStep()
            self.conv = True if iter >= self.__max_iter else False
            DictLogger("Convergence info",self.GetInfo())
        return self.conv

    def GetInfo(self) -> dict:
        info = {
            "type": "max_iter",
            "iter": f"{self.__optimization_problem.GetStep()} of {self.__max_iter}",
            "status": str("converged" if self.conv else "not converged")
        }
        return info

class L2ConvCriterion:
    @classmethod
    def GetDefaultParameters(cls):
        return Kratos.Parameters("""{
            "type"              : "l2_norm",
            "max_iter"          : 0,
            "tolerance"         : 1e-9
        }""")

    def __init__(self, parameters: Kratos.Parameters, optimization_problem: OptimizationProblem):
        parameters.ValidateAndAssignDefaults(self.GetDefaultParameters())
        self.__max_iter = parameters["max_iter"].GetInt()
        self.__optimization_problem = optimization_problem
        self.__tolerance = parameters["tolerance"].GetDouble()

    def IsConverged(self) -> bool:
        with TimeLogger("L2ConvCriterion::IsConverged", None, "Finished"):
            iter = self.__optimization_problem.GetStep()
            self.conv = True if iter >= self.__max_iter else False

            algorithm_buffered_data = ComponentDataView("algorithm", self.__optimization_problem).GetBufferedData()
            if not algorithm_buffered_data.HasValue("search_direction"):
                raise RuntimeError(f"Algorithm data does not contain computed \"search_direction\".\nData:\n{algorithm_buffered_data}" )

            self.norm = KratosOA.ExpressionUtils.NormL2(algorithm_buffered_data["search_direction"])
            if not self.conv:
                self.conv = True if self.norm <= self.__tolerance else False

            DictLogger("Convergence info",self.GetInfo())

        return self.conv

    def GetInfo(self) -> dict:
        info = {'type': 'l2_norm',
                'l2_norm': self.norm,
                'tolerance': self.__tolerance,
                'status': str("converged" if self.conv else "not converged")}
        return info
