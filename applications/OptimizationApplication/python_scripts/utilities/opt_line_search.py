import KratosMultiphysics as Kratos
import KratosMultiphysics.OptimizationApplication as KratosOA
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem import OptimizationProblem
from KratosMultiphysics.OptimizationApplication.utilities.component_data_view import ComponentDataView
from KratosMultiphysics.OptimizationApplication.utilities.logger_utilities import DictLogger
from KratosMultiphysics.OptimizationApplication.utilities.logger_utilities import time_decorator
from KratosMultiphysics.OptimizationApplication.algorithms.standardized_objective import StandardizedObjective
import math
import numpy as np


def CreateLineSearch(parameters: Kratos.Parameters, optimization_problem: OptimizationProblem, objective: StandardizedObjective):
    type = parameters["type"].GetString()
    if type == "const_step":
        return ConstStep(parameters, optimization_problem)
    elif type == "poly_step":
        return StepByPolynomialApproximation(parameters, optimization_problem, objective)
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
            raise RuntimeError(f"Algorithm data does not contain computed \"search_direction\".\nData:\n{algorithm_buffered_data}")

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
            self.step = self._init_step

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
            dy = KratosOA.ExpressionUtils.InnerProduct(d, y)
            dd = KratosOA.ExpressionUtils.InnerProduct(d, d)
            if not math.isclose(dy, 0.0, abs_tol=1e-16):
                self.unscaled_step = abs(dd / dy)
            else:
                self.unscaled_step = self._max_step

        if not math.isclose(norm, 0.0, abs_tol=1e-16):
            self.step = self.unscaled_step / norm
        else:
            self.step = self.unscaled_step

        if self.step > self._max_step:
            self.step = self._max_step

        DictLogger("Line Search info", self.GetInfo())

        return self.step

    def GetInfo(self) -> dict:
        info = {'type': 'BB_step',
                'unscaled_step': self.unscaled_step,
                'scaled_step': self.step,
                'max_step': self._max_step,
                'init_step': self._init_step}

        return info


class QNBBStep(BBStep):
    @time_decorator()
    def ComputeStep(self) -> KratosOA.CollectiveExpression:
        algorithm_buffered_data = ComponentDataView("algorithm", self._optimization_problem).GetBufferedData()
        norm = self.ComputeScaleFactor()
        self.step = algorithm_buffered_data.GetValue("search_direction", 0).Clone()
        self.step *= 0.0
        if not algorithm_buffered_data.HasValue("step_size"):
            algorithm_buffered_data["step_size"] = self.step
        self.step_numpy = self.step.Evaluate()
        if self._optimization_problem.GetStep() == 0:
            self.step_numpy[:] = self._init_step
        else:
            current_search_direction = algorithm_buffered_data.GetValue("search_direction", 0)
            previous_search_direction = algorithm_buffered_data.GetValue("search_direction", 1)
            y = previous_search_direction - current_search_direction
            y = y.Evaluate()
            d = algorithm_buffered_data.GetValue("control_field_update", 1)
            d = d.Evaluate()
            for i in range(len(y)):
                yy = y[i] * y[i]
                yd = y[i] * d[i]
                if math.isclose(yy, 0.0, abs_tol=1e-16):
                    self.step_numpy[i] = self._max_step
                else:
                    self.step_numpy[i] = abs(yd / yy)
                if self.step_numpy[i] > self._max_step:
                    self.step_numpy[i] = self._max_step

        if not math.isclose(norm, 0.0, abs_tol=1e-16):
            self.step_numpy[:] /= norm

        DictLogger("Line Search info", self.GetInfo())

        shape = [c.GetItemShape() for c in self.step.GetContainerExpressions()]
        KratosOA.CollectiveExpressionIO.Read(self.step, self.step_numpy, shape)
        return self.step

    def GetInfo(self) -> dict:
        info = {'type': 'QNBB_step',
                # 'max scaled_step': self.step.max(),
                'max_step': self._max_step,
                'init_step': self._init_step}

        return info


class StepByPolynomialApproximation(object):
    @classmethod
    def GetDefaultParameters(cls):
        return Kratos.Parameters("""{
            "type"                  : "poly_step",
            "init_step"             : 0,
            "gradient_scaling"      : "none",
            "degree_of_polynomial"  : 2,
            "step_damping"          : 1.0
        }""")

    def __init__(self, parameters: Kratos.Parameters, optimization_problem: OptimizationProblem, objective: StandardizedObjective):
        parameters.ValidateAndAssignDefaults(self.GetDefaultParameters())
        self.init_step = parameters["init_step"].GetDouble()
        self.degree_of_polynomial = parameters["degree_of_polynomial"].GetInt()
        self.step_damping = parameters["step_damping"].GetDouble()
        self.__optimization_problem = optimization_problem
        self.__gradient_scaling = parameters["gradient_scaling"].GetString()
        self.__objective = objective

    @time_decorator()
    def ComputeStep(self) -> KratosOA.CollectiveExpression:
        self.algorithm_buffered_data = ComponentDataView("algorithm", self.__optimization_problem).GetBufferedData()

        if not self.algorithm_buffered_data.HasValue("search_direction"):
            raise RuntimeError(f"Algorithm data does not contain computed \"search_direction\".\nData:\n{self.algorithm_buffered_data}")

        if self.__gradient_scaling == "inf_norm":
            norm = KratosOA.ExpressionUtils.NormInf(self.algorithm_buffered_data["search_direction"])
        elif self.__gradient_scaling == "l2_norm":
            norm = KratosOA.ExpressionUtils.NormL2(self.algorithm_buffered_data["search_direction"])
        elif self.__gradient_scaling == "none":
            norm = 1.0
        else:
            raise RuntimeError("\"gradient_scaling\" has unknown type.")

        if not math.isclose(norm, 0.0, abs_tol=1e-16):
            self.step_size = self.init_step / norm
        else:
            self.step_size = self.init_step

        self.step_size = self.find_best_step_size_from_polynomial_approximation(tests_step_size=self.init_step, degree_of_polynomial=self.degree_of_polynomial)

        DictLogger("Line Search info", self.GetInfo())

        return self.step_size

    def GetInfo(self) -> dict:
        info = {'type': 'constant',
                'unscaled_step': self.init_step,
                'scaled_step': self.step_size,
                'degree_polynomial': self.degree_of_polynomial}

        return info

    def find_best_step_size_from_polynomial_approximation(self, tests_step_size, degree_of_polynomial: int = 2):

        n__steps = degree_of_polynomial
        objective_values = np.zeros(n__steps+1)
        step_values = np.zeros(n__steps+1)
        update = self.algorithm_buffered_data["search_direction"] * tests_step_size
        self.__control_field = self.__objective.GetMasterControl().GetControlField()

        objective_values[0] = self.__objective.CalculateStandardizedValue(self.__control_field)
        step_values[0] = 0

        for i in range(1, n__steps+1):
            self.__control_field += update

            objective_values[i] = self.__objective.CalculateStandardizedValue(self.__control_field)
            step_values[i] = i * tests_step_size

        # Get coefficients of the approximation polynome
        polynomial_coefficients = np.polynomial.polynomial.polyfit(x=step_values, y=objective_values, deg=n__steps)

        # Calculate coefficients of the gradient polynome
        gradient_polynomial_coefficients = np.copy(polynomial_coefficients)
        for i in range(len(polynomial_coefficients)):
            gradient_polynomial_coefficients[i] *= i
        gradient_polynomial_coefficients = gradient_polynomial_coefficients[1:]

        if gradient_polynomial_coefficients.size == 1:
            best_step_length = -polynomial_coefficients[0]/polynomial_coefficients[1]
        elif gradient_polynomial_coefficients.size == 2:
            best_step_length = float(-(gradient_polynomial_coefficients[0] / (2*gradient_polynomial_coefficients[1])))
        else:
            roots = np.polynomial.polynomial.polyroots(polynomial_coefficients[1:])
            best_step_length = float(np.min(np.abs(roots)))

        if best_step_length > 16 * n__steps * tests_step_size:
            best_step_length = 16 * n__steps * tests_step_size

        best_step_length *= self.step_damping

        return best_step_length
