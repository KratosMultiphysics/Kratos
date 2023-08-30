import KratosMultiphysics as Kratos
import KratosMultiphysics.OptimizationApplication as KratosOA
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem import OptimizationProblem
from KratosMultiphysics.OptimizationApplication.utilities.component_data_view import ComponentDataView
<<<<<<<<< Temporary merge branch 1
=========
from KratosMultiphysics.OptimizationApplication.utilities.logger_utilities import DictLogger
from KratosMultiphysics.OptimizationApplication.utilities.logger_utilities import TimeLogger
import math
>>>>>>>>> Temporary merge branch 2

def CreateLineSearch(parameters: Kratos.Parameters, optimization_problem: OptimizationProblem, objective: StandardizedObjective):
    type = parameters["type"].GetString()
    if type == "const_step":
        return ConstStep(parameters, optimization_problem)
    else:
        raise RuntimeError(f"CreateConvergenceCriteria: unsupported convergence type {type}.")

class ConstStep(object):
    @classmethod
    def GetDefaultParameters(cls):
        return Kratos.Parameters("""{
            "type"               : "const_step",
            "init_step"          : 0,
            "gradient_scaling": ": inf_norm"
        }""")

    def __init__(self, parameters: Kratos.Parameters, optimization_problem: OptimizationProblem):
        parameters.ValidateAndAssignDefaults(self.GetDefaultParameters())
        self._init_step = parameters["init_step"].GetDouble()
        self._optimization_problem = optimization_problem
        self._gradient_scaling = parameters["gradient_scaling"].GetString()
        self.unscaled_step = self._init_step

    def ComputeStep(self) -> float:
<<<<<<<<< Temporary merge branch 1
        algorithm_buffered_data = ComponentDataView("algorithm", self.__optimization_problem).GetBufferedData()

        if not algorithm_buffered_data.HasValue("search_direction"):
            raise RuntimeError(f"Algorithm data does not contain computed \"search_direction\".\nData:\n{algorithm_buffered_data}" )

        if self.__gradient_scaling == "inf_norm":
            norm = KratosOA.ContainerExpressionUtils.NormInf(algorithm_buffered_data["search_direction"])
        elif self.__gradient_scaling == "l2_norm":
            norm = KratosOA.ContainerExpressionUtils.NormL2(algorithm_buffered_data["search_direction"])
        elif self.__gradient_scaling == "none":
            norm = 0.0
        else:
            raise RuntimeError("\"gradient_scaling\" has unknown type.")
        if norm:
            step = self.init_step / norm
        else:
            step =  self.init_step
        msg = f"""\t Line Search info:
            type          : constant
            unscaled_step : {self.init_step:0.6e}
            scaled_step   : {step:0.6e}"""
        print(msg)
        return step
=========
        with TimeLogger("ConstStep::ComputeStep", None, "Finished"):
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
            DictLogger("Line Search info", self.GetInfo())

        return self.step

    def GetInfo(self) -> dict:
        info = {'type': 'constant',
                'unscaled_step': self.unscaled_step,
                'scaled_step': self.step}
        return info

class BBStep(ConstStep):
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
        self._init_step = parameters["init_step"].GetDouble()
        self._max_step = parameters["max_step"].GetDouble()
        self._optimization_problem = optimization_problem
        self._gradient_scaling = parameters["gradient_scaling"].GetString()

    @time_decorator()
    def ComputeStep(self) -> KratosOA.CollectiveExpression:
        algorithm_buffered_data = ComponentDataView("algorithm", self._optimization_problem).GetBufferedData()
        norm = self.ComputeScaleFactor()
        if self._optimization_problem.GetStep() == 0:
            self.unscaled_step = self._init_step
        else:
            current_search_direction = algorithm_buffered_data.GetValue("search_direction", 0)
            previous_search_direction = algorithm_buffered_data.GetValue("search_direction", 1)
            y = previous_search_direction - current_search_direction
            d = algorithm_buffered_data.GetValue("control_field_update", 1)
            dy = KratosOA.ExpressionUtils.InnerProduct(d,y)
            dd = KratosOA.ExpressionUtils.InnerProduct(d,d)
            if not math.isclose(dy, 0.0, abs_tol=1e-16):
                self.unscaled_step = abs( dd / dy )
            else:
                self.unscaled_step = self._max_step

        if not math.isclose(norm, 0.0, abs_tol=1e-16):
            self.step = self.unscaled_step / norm
        else:
            self.step =  self.unscaled_step

        if self.step > self._max_step:
            self.step = self._max_step

        DictLogger("Line Search info",self.GetInfo())

        return self.step

    def GetInfo(self) -> dict:
        info = {'type': 'BB_step',
                'unscaled_step': self.unscaled_step,
                'scaled_step': self.step,
                'max_step': self._max_step,
                'init_step': self._init_step}

        return info
>>>>>>>>> Temporary merge branch 2

        return info