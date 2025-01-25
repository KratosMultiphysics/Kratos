import KratosMultiphysics as Kratos
import KratosMultiphysics.OptimizationApplication as KratosOA
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem import OptimizationProblem
from KratosMultiphysics.OptimizationApplication.utilities.component_data_view import ComponentDataView
from KratosMultiphysics.OptimizationApplication.utilities.logger_utilities import DictLogger
from KratosMultiphysics.OptimizationApplication.utilities.logger_utilities import time_decorator
import numpy

import math

def CreateLineSearch(parameters: Kratos.Parameters, optimization_problem: OptimizationProblem):
    type = parameters["type"].GetString()
    if type == "const_step":
        return ConstStep(parameters, optimization_problem)
    elif type == "BB_step":
        return BBStep(parameters, optimization_problem)
    elif type == "QNBB_step":
        return QNBBStep(parameters, optimization_problem)
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
        self._init_step = parameters["init_step"].GetDouble()
        self._optimization_problem = optimization_problem
        self._gradient_scaling = parameters["gradient_scaling"].GetString()
        self.unscaled_step = self._init_step

    def ComputeScaleFactor(self) -> float:
        algorithm_buffered_data = ComponentDataView("algorithm", self._optimization_problem).GetBufferedData()

        if not algorithm_buffered_data.HasValue("search_direction"):
            raise RuntimeError(f"Algorithm data does not contain computed \"search_direction\".\nData:\n{algorithm_buffered_data}" )

        if self._gradient_scaling == "inf_norm":
            norm = KratosOA.ExpressionUtils.NormInf(algorithm_buffered_data["search_direction"])
        elif self._gradient_scaling == "l2_norm":
            norm = KratosOA.ExpressionUtils.NormL2(algorithm_buffered_data["search_direction"])
        elif self._gradient_scaling == "none":
            norm = 1.0
        else:
            raise RuntimeError("\"gradient_scaling\" has unknown type.")

        return norm

    @time_decorator()
    def ComputeStep(self) -> float:
        norm = self.ComputeScaleFactor()
        if not math.isclose(norm, 0.0, abs_tol=1e-16):
            self.step = self._init_step / norm
        else:
            self.step =  self._init_step

        DictLogger("Line Search info",self.GetInfo())

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
        if math.isclose(norm, 0.0, abs_tol=1e-16):
            norm = 1.0
        if self._optimization_problem.GetStep() == 0:
            self.step = self._init_step / norm
        else:
            current_search_direction = algorithm_buffered_data.GetValue("search_direction", 0)
            previous_search_direction = algorithm_buffered_data.GetValue("search_direction", 1)
            y = previous_search_direction - current_search_direction
            d = algorithm_buffered_data.GetValue("control_field_update", 1)
            dy = KratosOA.ExpressionUtils.InnerProduct(d,y)
            dd = KratosOA.ExpressionUtils.InnerProduct(d,d)
            if math.isclose(dy, 0.0, abs_tol=1e-16):
                self.step = self._max_step / norm
            else:
                self.step = abs( dd / dy )

            if self.step > self._max_step / norm:
                self.step = self._max_step / norm

        DictLogger("Line Search info",self.GetInfo())

        return self.step

    def GetInfo(self) -> dict:
        info = {'type': 'BB_step',
                'step': self.step,
                'max_step': self._max_step,
                'init_step': self._init_step}

        return info

class QNBBStep(BBStep):
    @time_decorator()
    def ComputeStep(self) -> KratosOA.CollectiveExpression:
        algorithm_buffered_data = ComponentDataView("algorithm", self._optimization_problem).GetBufferedData()
        norm = self.ComputeScaleFactor()
        if math.isclose(norm, 0.0, abs_tol=1e-16):
            norm = 1.0
        self.step = algorithm_buffered_data.GetValue("search_direction", 0).Clone()
        self.step *= 0.0
        if not algorithm_buffered_data.HasValue("step_size"):
            algorithm_buffered_data["step_size"] = self.step
        self.step_numpy = self.step.Evaluate()
        if self._optimization_problem.GetStep() == 0:
            self.step_numpy[:] = self._init_step / norm
        else:
            current_search_direction = algorithm_buffered_data.GetValue("search_direction", 0)
            previous_search_direction = algorithm_buffered_data.GetValue("search_direction", 1)
            y = previous_search_direction - current_search_direction
            y = y.Evaluate()
            d = algorithm_buffered_data.GetValue("control_field_update", 1)
            d = d.Evaluate()
            for i in range(len(y)):
                dy = numpy.dot(d[i], y[i])
                dd = numpy.dot(d[i], d[i])
                
                if math.isclose(dy, 0.0, abs_tol=1e-16):
                    self.step_numpy[i] = self._max_step / norm
                else:
                    self.step_numpy[i] = abs( dd / dy )

                if isinstance(d[i], (float, int, numpy.float64)) and self.step_numpy[i] > self._max_step / norm:
                    self.step_numpy[i] = self._max_step / norm
                elif isinstance(d[i], (numpy.ndarray)) and self.step_numpy[i][0] > self._max_step / norm:
                        self.step_numpy[i] = self._max_step / norm
                
        DictLogger("Line Search info",self.GetInfo())

        shape = [c.GetItemShape() for c in self.step.GetContainerExpressions()]
        KratosOA.CollectiveExpressionIO.Read(self.step, self.step_numpy, shape)
        return self.step

    def GetInfo(self) -> dict:
        info = {'type': 'QNBB_step',
                # 'max scaled_step': self.step.max(),
                'max_step': self._max_step,
                'init_step': self._init_step}

        return info