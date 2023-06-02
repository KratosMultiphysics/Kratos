import KratosMultiphysics as Kratos
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem import OptimizationProblem
from KratosMultiphysics.OptimizationApplication.utilities.component_data_view import ComponentDataView
import KratosMultiphysics.OptimizationApplication as KratosOA

def CreateConvergenceCriteria(parameters: Kratos.Parameters, optimization_problem: OptimizationProblem):
    type = parameters["type"].GetString()
    if type == "max_iter":
        return MaxIterConvCriterium(parameters, optimization_problem)
    elif type == "l2_norm":
        return L2ConvCriterium(parameters, optimization_problem)
    else:
        raise RuntimeError(f"CreateConvergenceCriteria: unsupported convergence type {type}.")

class MaxIterConvCriterium:
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
        iter = self.__optimization_problem.GetStep()
        conv = True if iter >= self.__max_iter else False
        msg = f"""\t Convergence info: 
            type          : max_iter
            iter          : {iter} of {self.__max_iter}
            status        : {"converged" if conv else "not converged"}"""
        print(msg)
        return conv
    
class L2ConvCriterium:
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
        iter = self.__optimization_problem.GetStep()
        conv = True if iter >= self.__max_iter else False

        algorithm_buffered_data = ComponentDataView("algorithm", self.__optimization_problem).GetBufferedData()
        if not algorithm_buffered_data.HasValue("search_direction"):
            raise RuntimeError(f"Algorithm data does not contain computed \"search_direction\".\nData:\n{algorithm_buffered_data}" )

        norm = KratosOA.ContainerExpressionUtils.NormL2(algorithm_buffered_data["search_direction"])
        if not conv:
            conv = True if norm <= self.__tolerance else False 

        msg = f"""\t Convergence info: 
            type          : l2_norm
            l2_norm       : {norm:0.6e}
            tolerance     : {self.__tolerance:0.6e}
            iter          : {iter} of {self.__max_iter}
            status        : {"converged" if conv else "not converged"}"""
        print(msg)
        return conv
