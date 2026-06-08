import KratosMultiphysics as Kratos
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

        search_direction: Kratos.TensorAdaptors.DoubleCombinedTensorAdaptor = algorithm_buffered_data["search_direction"]
        if self._gradient_scaling == "inf_norm":
            return numpy.max(numpy.abs(search_direction.data))
        elif self._gradient_scaling == "l2_norm":
            return numpy.linalg.norm(search_direction.data)
        elif self._gradient_scaling == "none":
            return 1.0
        else:
            raise RuntimeError("\"gradient_scaling\" has unknown type.")

    @time_decorator()
    def ComputeStep(self) -> Kratos.TensorAdaptors.DoubleCombinedTensorAdaptor:
        norm = self.ComputeScaleFactor()
        if not math.isclose(norm, 0.0, abs_tol=1e-16):
            self.step = self._init_step / norm
        else:
            self.step =  self._init_step

        DictLogger("Line Search info",self.GetInfo())

        algorithm_buffered_data = ComponentDataView("algorithm", self._optimization_problem).GetBufferedData()
        step = Kratos.TensorAdaptors.DoubleCombinedTensorAdaptor(algorithm_buffered_data["search_direction"])
        step.data[:] = self.step
        return step

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
    def ComputeStep(self) -> Kratos.TensorAdaptors.DoubleCombinedTensorAdaptor:
        algorithm_buffered_data = ComponentDataView("algorithm", self._optimization_problem).GetBufferedData()
        norm = self.ComputeScaleFactor()
        if math.isclose(norm, 0.0, abs_tol=1e-16):
            norm = 1.0
        if self._optimization_problem.GetStep() == 0:
            self.step = self._init_step / norm
        else:
            current_search_direction: Kratos.TensorAdaptors.DoubleCombinedTensorAdaptor = algorithm_buffered_data.GetValue("search_direction", 0)
            previous_search_direction: Kratos.TensorAdaptors.DoubleCombinedTensorAdaptor = algorithm_buffered_data.GetValue("search_direction", 1)
            y = previous_search_direction.data - current_search_direction.data
            d: Kratos.TensorAdaptors.DoubleCombinedTensorAdaptor = algorithm_buffered_data.GetValue("control_field_update", 1)
            dy = numpy.inner(d.data[:].ravel(), y.ravel())
            dd = numpy.inner(d.data[:].ravel(), d.data[:].ravel())
            if math.isclose(dy, 0.0, abs_tol=1e-16):
                self.step = self._max_step / norm
            else:
                self.step = abs( dd / dy )

            if self.step > self._max_step / norm:
                self.step = self._max_step / norm

        DictLogger("Line Search info",self.GetInfo())

        algorithm_buffered_data = ComponentDataView("algorithm", self._optimization_problem).GetBufferedData()
        step = Kratos.TensorAdaptors.DoubleCombinedTensorAdaptor(algorithm_buffered_data["search_direction"])
        step.data[:] = self.step
        return step

    def GetInfo(self) -> dict:
        info = {'type': 'BB_step',
                'step': self.step,
                'max_step': self._max_step,
                'init_step': self._init_step}

        return info

class QNBBStep(BBStep):
    @time_decorator()
    def ComputeStep(self) -> Kratos.TensorAdaptors.DoubleCombinedTensorAdaptor:
        algorithm_buffered_data = ComponentDataView("algorithm", self._optimization_problem).GetBufferedData()
        norm = self.ComputeScaleFactor()
        if math.isclose(norm, 0.0, abs_tol=1e-16):
            norm = 1.0
        self.step = Kratos.TensorAdaptors.DoubleCombinedTensorAdaptor(algorithm_buffered_data.GetValue("search_direction", 0), perform_store_data_recursively=False)
        self.step.data[:] = 0.0
        if not algorithm_buffered_data.HasValue("step_size"):
            algorithm_buffered_data["step_size"] = self.step

        if self._optimization_problem.GetStep() == 0:
            self.step.data[:] = self._init_step / norm
        else:
            current_search_direction: Kratos.TensorAdaptors.DoubleCombinedTensorAdaptor = algorithm_buffered_data.GetValue("search_direction", 0)
            previous_search_direction: Kratos.TensorAdaptors.DoubleCombinedTensorAdaptor = algorithm_buffered_data.GetValue("search_direction", 1)
            y = previous_search_direction.data - current_search_direction.data
            d: Kratos.TensorAdaptors.DoubleCombinedTensorAdaptor = algorithm_buffered_data.GetValue("control_field_update", 1)
            d = d.data[:]
            for i in range(len(y)):
                dy = numpy.dot(d[i], y[i])
                dd = numpy.dot(d[i], d[i])

                if math.isclose(dy, 0.0, abs_tol=1e-16):
                    self.step.data[i] = self._max_step / norm
                else:
                    self.step.data[i] = abs( dd / dy )

                if isinstance(d[i], (float, int, numpy.float64)) and self.step.data[i] > self._max_step / norm:
                    self.step.data[i] = self._max_step / norm
                elif isinstance(d[i], (numpy.ndarray)) and self.step.data[i][0] > self._max_step / norm:
                        self.step.data[i] = self._max_step / norm

        DictLogger("Line Search info",self.GetInfo())
        self.step.StoreData()
        return self.step

    def GetInfo(self) -> dict:
        info = {'type': 'QNBB_step',
                # 'max scaled_step': self.step.max(),
                'max_step': self._max_step,
                'init_step': self._init_step}

        return info