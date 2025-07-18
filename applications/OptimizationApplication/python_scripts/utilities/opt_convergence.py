import abc, numpy
import KratosMultiphysics as Kratos
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem import OptimizationProblem
from KratosMultiphysics.OptimizationApplication.utilities.component_data_view import ComponentDataView
import KratosMultiphysics.OptimizationApplication as KratosOA
from KratosMultiphysics.OptimizationApplication.utilities.logger_utilities import DictLogger
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
    elif type == "magnitude_reduction":
        return MagnitudeReductionCriterion(parameters, optimization_problem)
    elif type == "plateau":
        return PlateauConvergenceCriteria(parameters, optimization_problem)
    else:
        raise RuntimeError(f"CreateConvergenceCriteria: unsupported convergence type {type}.")

class ConvergenceCriterion(abc.ABC):
    """Base class for the convergence criteria.

    This is the base class convergence criteria which is used to determine
    whether the objective has converged or not.

    This convergence criteria always require max number of iterations to
    be specified to break the optimization no matter the convergence is reached
    or not.

    """
    def __init__(self, max_iterations: int, optimization_problem: OptimizationProblem) -> None:
        self.__max_iterations = max_iterations
        self.__optimization_problem = optimization_problem

    @abc.abstractmethod
    def IsConverged(self) -> bool:
        """Returns whether the objective is converged or not."""
        pass

    def GetMaxIterations(self) -> int:
        return self.__max_iterations

    def IsMaxIterationsReached(self) -> bool:
        return self.__optimization_problem.GetStep() >= self.GetMaxIterations()

class MaxIterConvCriterion(ConvergenceCriterion):
    """
    A convergence criterion class that checks if the maximum number of iterations has been reached.

    Methods
    -------
    GetDefaultParameters():
        Returns the default parameters for the convergence criterion.

    __init__(parameters: Kratos.Parameters, optimization_problem: OptimizationProblem):
        Initializes the convergence criterion with the given parameters and optimization problem.

    IsConverged(search_direction=None) -> bool:
        Checks if the optimization problem has converged based on the maximum number of iterations.

    GetInfo() -> dict:
        Returns a dictionary with information about the current state of convergence.
    """
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
        super().__init__(self.__max_iter, self.__optimization_problem)

    @time_decorator()
    def IsConverged(self) -> bool:
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

class L2ConvCriterion(ConvergenceCriterion):
    """
    A class to check the convergence of an optimization problem using the L2 norm criterion.

    Methods
    -------
    GetDefaultParameters():
        Returns the default parameters for the L2ConvCriterion.

    __init__(parameters: Kratos.Parameters, optimization_problem: OptimizationProblem):
        Initializes the L2ConvCriterion with given parameters and optimization problem.

    IsConverged() -> bool:
        Checks if the optimization problem has converged based on the L2 norm of the search direction.

    GetInfo() -> dict:
        Returns a dictionary with information about the current state of convergence.
    """
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
        super().__init__(self.__max_iter, self.__optimization_problem)

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

class AverageAbsoluteImprovement(ConvergenceCriterion):
    """
    A class to monitor the convergence of an optimization problem based on the average absolute improvement
    of a tracked value over a specified number of iterations.

    Methods
    -------
    GetDefaultParameters():
        Returns the default parameters for the convergence criteria.

    __init__(parameters: Kratos.Parameters, optimization_problem: OptimizationProblem):
        Initializes the convergence criteria with the given parameters and optimization problem.

    IsConverged() -> bool:
        Checks if the optimization problem has converged based on the average absolute improvement.

    GetInfo() -> dict:
        Returns a dictionary with information about the current state of the convergence criteria.
    """
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
        super().__init__(self.__max_iter, self.__optimization_problem)

    def IsMaxIterationsReached(self):
        return self.__optimization_problem.GetStep() >= self.__max_iter

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

class TargetValueCriterion(ConvergenceCriterion):
    """
    A class to check the convergence of an optimization problem based on a target value criterion.

    Methods
    -------
    GetDefaultParameters():
        Returns the default parameters for the criterion.

    __init__(parameters: Kratos.Parameters, optimization_problem: OptimizationProblem):
        Initializes the criterion with given parameters and optimization problem.

    IsConverged() -> bool:
        Checks if the optimization problem has converged based on the target value.

    GetInfo() -> dict:
        Returns a dictionary with information about the current state of convergence.
    """
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
        super().__init__(self.__max_iter, self.__optimization_problem)

    def IsMaxIterationsReached(self):
        return self.__optimization_problem.GetStep() >= self.__max_iter

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

class MagnitudeReductionCriterion(ConvergenceCriterion):
    @classmethod
    def GetDefaultParameters(cls):
        return Kratos.Parameters("""{
            "type"                            : "magnitude_reduction",
            "max_iter"                        : 0,
            "target_scaling_factor"           : 1e-3,
            "machine_precision"               : 1e-16
        }""")

    def __init__(self, parameters: Kratos.Parameters, optimization_problem: OptimizationProblem):
        parameters.ValidateAndAssignDefaults(self.GetDefaultParameters())
        self.__max_iter = parameters["max_iter"].GetInt()
        self.__optimization_problem = optimization_problem
        self.__target_scaling_factor = parameters["target_scaling_factor"].GetDouble()
        self.__machine_precision = parameters["machine_precision"].GetDouble()
        super().__init__(self.__max_iter, self.__optimization_problem)

    @time_decorator()
    def IsConverged(self, ) -> bool:
        iter = self.__optimization_problem.GetStep()
        self.conv = iter >= self.__max_iter

        algorithm_buffered_data = ComponentDataView("algorithm", self.__optimization_problem).GetBufferedData()
        if not algorithm_buffered_data.HasValue("std_obj_value"):
            raise RuntimeError(f"Algorithm data does not contain computed \"std_obj_value\".\nData:\n{algorithm_buffered_data}" )
        self.value = algorithm_buffered_data["std_obj_value"]

        if iter == 0:
            self.__target_value = max(self.value * self.__target_scaling_factor, self.__machine_precision)
        if not self.conv:
            self.conv = self.value <= self.__target_value

        DictLogger("Convergence info",self.GetInfo())

        return self.conv

    def GetInfo(self) -> dict:
        info = {'type': 'magnitude_reduction',
                'current_value': self.value,
                'target_value': self.__target_value,
                "iter": f"{self.__optimization_problem.GetStep()} of {self.__max_iter}",
                'status': str("converged" if self.conv else "not converged")}
        return info

class PlateauConvergenceCriteria(ConvergenceCriterion):
    @classmethod
    def GetDefaultParameters(cls):
        return Kratos.Parameters("""{
            "type"             : "plateau",
            "max_iter"         : 500,
            "min_iter"         : 100,
            "min_delta"        : 1e-8,
            "patience"         : 50,
            "is_delta_relative": true
        }""")

    def __init__(self, parameters: Kratos.Parameters, optimization_problem: OptimizationProblem):
        parameters.ValidateAndAssignDefaults(self.GetDefaultParameters())
        super().__init__(parameters["max_iter"].GetInt(), optimization_problem)

        self.optimization_problem = optimization_problem
        self.min_delta = parameters["min_delta"].GetDouble()
        self.min_iter = parameters["min_iter"].GetInt()

        self.history = Kratos.Vector(parameters["patience"].GetInt(), 0.0)
        if self.history.Size() == 0:
            raise RuntimeError(f"The patience should be a positive integer [ patience = {self.history.Size()} ].")

        self.cyclic_index = 0
        self.last_step = -1
        self.scaling_factor = None
        if not parameters["is_delta_relative"].GetBool():
            self.scaling_factor = 1.0

        self.values = []

    @time_decorator()
    def IsConverged(self):
        algorithm_buffered_data = ComponentDataView("algorithm", self.optimization_problem).GetBufferedData()
        if not algorithm_buffered_data.HasValue("std_obj_value"):
            raise RuntimeError(f"Algorithm data does not contain computed \"std_obj_value\".\nData:\n{algorithm_buffered_data}" )
        self.value = algorithm_buffered_data["std_obj_value"]

        current_step = self.optimization_problem.GetStep()
        if current_step > self.last_step:
            self.cyclic_index = (self.cyclic_index + 1) % self.history.Size()
            self.last_step = current_step

            if self.scaling_factor is None:
                self.scaling_factor = self.value

        self.history[self.cyclic_index] = self.value

        # self.values.append(self.value)
        # if self.optimization_problem.GetStep() > self.min_iter:
        #     n = len(self.values)
        #     fft_vals = numpy.fft.fft(y)
        #     fft_freqs = numpy.fft.fftfreq(n, d=1)  # d=1 since STEP increments by 1
        #     pos_mask = fft_freqs >= 0
        #     fft_freqs = fft_freqs[pos_mask]
        #     fft_magnitude = np.abs(fft_vals[pos_mask])
        # else:
        #     return False

        # now calculate the max and min in the history
        self.max_value = max(self.history)
        self.min_value = min(self.history)
        self.patience_max_delta = (self.max_value - self.min_value) / self.scaling_factor
        self.conv = self.patience_max_delta < self.min_delta

        DictLogger("Convergence info",self.GetInfo())

        return self.min_iter < self.optimization_problem.GetStep() and self.conv

    def GetInfo(self) -> dict:
        info = {'type': 'plateau',
                'current_value': self.value,
                'patience min_value': self.min_value,
                'patience max_value': self.max_value,
                'patience max_delta': self.patience_max_delta,
                'patience target max_delta': self.min_delta,
                'scaling factor'    : self.scaling_factor,
                'minimum number of iters': self.min_iter,
                "iter": f"{self.optimization_problem.GetStep()} of {self.GetMaxIterations()}",
                'status': str("converged" if self.conv else "not converged")}
        return info