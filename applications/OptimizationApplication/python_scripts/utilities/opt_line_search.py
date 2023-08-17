import KratosMultiphysics as Kratos
import KratosMultiphysics.OptimizationApplication as KratosOA
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem import OptimizationProblem
from KratosMultiphysics.OptimizationApplication.utilities.component_data_view import ComponentDataView
from KratosMultiphysics.OptimizationApplication.utilities.logger_utilities import DictLogger
from KratosMultiphysics.OptimizationApplication.utilities.logger_utilities import time_decorator

import math

def CreateLineSearch(parameters: Kratos.Parameters, optimization_problem: OptimizationProblem):
    type = parameters["type"].GetString()
    if type == "const_step":
        return ConstStep(parameters, optimization_problem)
    elif type == "BB_step":
        return BBStep(parameters, optimization_problem)
    else:
        raise RuntimeError(f"CreateConvergenceCriteria: unsupported convergence type {type}.")

class ConstStep():
    @classmethod
    def GetDefaultParameters(cls):
        return Kratos.Parameters("""{
            "type"               : "const_step",
            "init_step"          : 0,
            "gradient_scaling": ": inf_norm"
        }""")

    def __init__(self, parameters: Kratos.Parameters, optimization_problem: OptimizationProblem):
        parameters.ValidateAndAssignDefaults(self.GetDefaultParameters())
        self.init_step = parameters["init_step"].GetDouble()
        self.__optimization_problem = optimization_problem
        self.__gradient_scaling = parameters["gradient_scaling"].GetString()

    @time_decorator()
    def ComputeStep(self) -> float:
        algorithm_buffered_data = ComponentDataView("algorithm", self.__optimization_problem).GetBufferedData()

        if not algorithm_buffered_data.HasValue("search_direction"):
            raise RuntimeError(f"Algorithm data does not contain computed \"search_direction\".\nData:\n{algorithm_buffered_data}" )

        if self.__gradient_scaling == "inf_norm":
            norm = KratosOA.ExpressionUtils.NormInf(algorithm_buffered_data["search_direction"])
        elif self.__gradient_scaling == "l2_norm":
            norm = KratosOA.ExpressionUtils.NormL2(algorithm_buffered_data["search_direction"])
        elif self.__gradient_scaling == "none":
            norm = 1.0
        else:
            raise RuntimeError("\"gradient_scaling\" has unknown type.")
        if not math.isclose(norm, 0.0, abs_tol=1e-16):
            self.step = self.init_step / norm
        else:
            self.step =  self.init_step

        DictLogger("Line Search info",self.GetInfo())

        return self.step

    def GetInfo(self) -> dict:
        info = {'type': 'constant',
                'unscaled_step': self.init_step,
                'scaled_step': self.step}

        return info
    
class BBStep():
    @classmethod
    def GetDefaultParameters(cls):
        return Kratos.Parameters("""{
            "type"               : "BB_step",
            "init_step"          : 0,
            "max_step"           : 0,
            "gradient_scaling"   : "inf_norm"
        }""")

    def __init__(self, parameters: Kratos.Parameters, optimization_problem: OptimizationProblem):
        parameters.ValidateAndAssignDefaults(self.GetDefaultParameters())
        self.init_step = parameters["init_step"].GetDouble()
        self.max_step = parameters["max_step"].GetDouble()
        self.__optimization_problem = optimization_problem
        self.__gradient_scaling = parameters["gradient_scaling"].GetString()

    @time_decorator()
    def ComputeStep(self) -> float:
        algorithm_buffered_data = ComponentDataView("algorithm", self.__optimization_problem).GetBufferedData()

        if not algorithm_buffered_data.HasValue("search_direction"):
            raise RuntimeError(f"Algorithm data does not contain computed \"search_direction\".\nData:\n{algorithm_buffered_data}" )

        if self.__gradient_scaling == "inf_norm":
            norm = KratosOA.ExpressionUtils.NormInf(algorithm_buffered_data["search_direction"])
        elif self.__gradient_scaling == "l2_norm":
            norm = KratosOA.ExpressionUtils.NormL2(algorithm_buffered_data["search_direction"])
        elif self.__gradient_scaling == "none":
            norm = 1.0
        else:
            raise RuntimeError("\"gradient_scaling\" has unknown type.")
        
        if self.__optimization_problem.GetStep() == 0:
            self.unscaled_step = self.init_step
        else:
            current_search_direction = algorithm_buffered_data.GetValue("search_direction", 0)
            previous_search_direction = algorithm_buffered_data.GetValue("search_direction", 1)
            y = previous_search_direction - current_search_direction
            d = algorithm_buffered_data.GetValue("control_field_update", 1)
            dy = KratosOA.ExpressionUtils.InnerProduct(d,y)
            dd = KratosOA.ExpressionUtils.InnerProduct(d,d)
            if not math.isclose(dy, 0.0, abs_tol=1e-16):
                self.unscaled_step = abs(dd / dy)
            if math.isclose(dy, 0.0, abs_tol=1e-16):
                self.unscaled_step = self.max_step
            if math.isclose(dd, 0.0, abs_tol=1e-16):
                self.unscaled_step = 0.0

        if not math.isclose(norm, 0.0, abs_tol=1e-16):
            self.step = self.unscaled_step / norm
        else:
            self.step =  self.unscaled_step

        if self.step > self.max_step:
            self.step = self.max_step

        DictLogger("Line Search info",self.GetInfo())

        return self.step

    def GetInfo(self) -> dict:
        info = {'type': 'BB_step',
                'unscaled_step': self.unscaled_step,
                'scaled_step': self.step,
                'max_step': self.max_step,
                'init_step': self.init_step}

        return info
