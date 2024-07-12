import KratosMultiphysics as Kratos
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem import OptimizationProblem
from KratosMultiphysics.OptimizationApplication.utilities.component_data_view import ComponentDataView
import KratosMultiphysics.OptimizationApplication as KratosOA
from KratosMultiphysics.OptimizationApplication.utilities.logger_utilities import DictLogger
from KratosMultiphysics.OptimizationApplication.utilities.logger_utilities import TimeLogger
from KratosMultiphysics.OptimizationApplication.utilities.logger_utilities import time_decorator

def CreateConvergenceCriteria(parameters: Kratos.Parameters, optimization_problem: OptimizationProblem):
    type = parameters["type"].GetString()
    if type == "max_iter":
        return MaxIterConvCriterion(parameters, optimization_problem)
    elif type == "l2_norm":
        return L2ConvCriterion(parameters, optimization_problem)
    elif type == "aver_abs_delta":
        return AverageAbsoluteImprovement(parameters, optimization_problem)
    elif type == "target_value":
        return TargetValueCriterion(parameters, optimization_problem)
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

    @time_decorator()
    def IsConverged(self, search_direction=None) -> bool:
        iter = self.__optimization_problem.GetStep()
        self.conv = iter >= self.__max_iter
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

    @time_decorator()
    def IsConverged(self) -> bool:
        iter = self.__optimization_problem.GetStep()
        self.conv = iter >= self.__max_iter

        algorithm_buffered_data = ComponentDataView("algorithm", self.__optimization_problem).GetBufferedData()
        if not algorithm_buffered_data.HasValue("search_direction"):
            raise RuntimeError(f"Algorithm data does not contain computed \"search_direction\".\nData:\n{algorithm_buffered_data}" )

        self.norm = KratosOA.ExpressionUtils.NormL2(algorithm_buffered_data["search_direction"])
        if not self.conv:
            self.conv = self.norm <= self.__tolerance

        DictLogger("Convergence info",self.GetInfo())

        return self.conv

    def GetInfo(self) -> dict:
        info = {'type': 'l2_norm',
                'l2_norm': self.norm,
                'tolerance': self.__tolerance,
                "iter": f"{self.__optimization_problem.GetStep()} of {self.__max_iter}",
                'status': str("converged" if self.conv else "not converged")}
        return info

class AverageAbsoluteImprovement:
    @classmethod
    def GetDefaultParameters(cls):
        return Kratos.Parameters("""{
            "type"              : "aver_abs_delta",
            "max_iter"          : 0,
            "tracked_iter"      : 5,
            "tolerance"         : 1e-6
        }""")

    def __init__(self, parameters: Kratos.Parameters, optimization_problem: OptimizationProblem):
        parameters.ValidateAndAssignDefaults(self.GetDefaultParameters())
        self.__max_iter = parameters["max_iter"].GetInt()
        self.__tracked_iter = parameters["tracked_iter"].GetInt()
        self.__optimization_problem = optimization_problem
        self.__tolerance = parameters["tolerance"].GetDouble()
        self.__abs_value_change = []
        self.value = None

    @time_decorator()
    def IsConverged(self) -> bool:
        iter = self.__optimization_problem.GetStep()
        self.conv = iter >= self.__max_iter

        step = self.__optimization_problem.GetStep()
        algorithm_buffered_data = ComponentDataView("algorithm", self.__optimization_problem).GetBufferedData()
        if step == 0:
            self.__init_value = algorithm_buffered_data.GetValue("std_obj_value", 0)
            self.conv = False
        if step > 0:
            if not algorithm_buffered_data.HasValue("std_obj_value"):
                raise RuntimeError(f"Algorithm data does not contain computed \"std_obj_value\".\nData:\n{algorithm_buffered_data}" )
            value = algorithm_buffered_data.GetValue("std_obj_value", 0)
            prev_value = algorithm_buffered_data.GetValue("std_obj_value", 1)
            self.__abs_value_change.append((value - prev_value)/self.__init_value)

        if step >= self.__tracked_iter - 1:
            self.value = 0.0
            for i in range(self.__tracked_iter - 1):
                self.value += self.__abs_value_change[-1 - i]
            self.value = abs(self.value) / (self.__tracked_iter - 1) 
            if not self.conv:
                self.conv = self.value <= self.__tolerance

        DictLogger("Convergence info",self.GetInfo())

        return self.conv

    def GetInfo(self) -> dict:
        info = {'type': 'aver_abs_delta',
                'aver_abs_delta': str(self.value),
                'tolerance': self.__tolerance,
                'tracked_iter': self.__tracked_iter,
                "iter": f"{self.__optimization_problem.GetStep()} of {self.__max_iter}",
                'status': str("converged" if self.conv else "not converged")}
        return info
    
class TargetValueCriterion:
    @classmethod
    def GetDefaultParameters(cls):
        return Kratos.Parameters("""{
            "type"              : "target_value",
            "max_iter"          : 0,
            "target_value"      : 1e-9
        }""")

    def __init__(self, parameters: Kratos.Parameters, optimization_problem: OptimizationProblem):
        parameters.ValidateAndAssignDefaults(self.GetDefaultParameters())
        self.__max_iter = parameters["max_iter"].GetInt()
        self.__optimization_problem = optimization_problem
        self.__target_value = parameters["target_value"].GetDouble()

    @time_decorator()
    def IsConverged(self) -> bool:
        iter = self.__optimization_problem.GetStep()
        self.conv = iter >= self.__max_iter

        algorithm_buffered_data = ComponentDataView("algorithm", self.__optimization_problem).GetBufferedData()
        if not algorithm_buffered_data.HasValue("std_obj_value"):
            raise RuntimeError(f"Algorithm data does not contain computed \"std_obj_value\".\nData:\n{algorithm_buffered_data}" )
        self.value = algorithm_buffered_data["std_obj_value"]
        if not self.conv:
            self.conv = self.value <= self.__target_value

        DictLogger("Convergence info",self.GetInfo())

        return self.conv

    def GetInfo(self) -> dict:
        info = {'type': 'value',
                'current_value': self.value,
                'target_value': self.__target_value,
                "iter": f"{self.__optimization_problem.GetStep()} of {self.__max_iter}",
                'status': str("converged" if self.conv else "not converged")}
        return info