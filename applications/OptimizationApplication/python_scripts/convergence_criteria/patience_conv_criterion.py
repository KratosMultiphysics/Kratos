import typing
import numpy
import KratosMultiphysics as Kratos
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem import OptimizationProblem
from KratosMultiphysics.OptimizationApplication.utilities.component_data_view import ComponentDataView
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem_utilities import GetComponentHavingDataByFullName
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem_utilities import GetComponentValueByFullName
from KratosMultiphysics.OptimizationApplication.convergence_criteria.convergence_criterion import ConvergenceCriterion
from KratosMultiphysics.OptimizationApplication.utilities.logger_utilities import time_decorator

def Factory(_: Kratos.Model, parameters: Kratos.Parameters, optimization_problem: OptimizationProblem) -> ConvergenceCriterion:
    if not parameters.Has("settings"):
        raise RuntimeError(f"PatienceConvCriterion instantiation requires a \"settings\" in parameters [ parameters = {parameters}].")
    return PatienceConvCriterion(parameters["settings"], optimization_problem)

class PatienceConvCriterion(ConvergenceCriterion):
    """
    Convergence criterion based on the 'patience' concept, commonly used in optimization and machine learning.

    This criterion monitors the improvement of a specified value (e.g., objective function) over iterations.
    Convergence is declared if no significant improvement (greater than a specified tolerance) is observed
    for a consecutive number of iterations (patience_itr), after a minimum number of iterations (minimum_itr)
    have been performed.

    Raises:
        RuntimeError: If tolerance is negative, minimum_itr is not positive, patience_itr is not greater than minimum_itr,
                      or the monitored value is not a float or int.
    """
    @classmethod
    def GetDefaultParameters(cls):
        return Kratos.Parameters("""{
            "component_name": "algorithm",
            "value_name"    : "std_obj_value",
            "minimum_itr"   : 100,
            "patience_itr"  : 50,
            "tolerance"     : 1e-6
        }""")

    def __init__(self, parameters: Kratos.Parameters, optimization_problem: OptimizationProblem):
        parameters.ValidateAndAssignDefaults(self.GetDefaultParameters())

        self.__component_name = parameters["component_name"].GetString()
        self.__value_name = parameters["value_name"].GetString()
        self.__tolerance = parameters["tolerance"].GetDouble()
        self.__patience_itr = parameters["patience_itr"].GetInt()
        self.__minimum_itr = parameters["minimum_itr"].GetInt()
        self.__optimization_problem = optimization_problem
        self.__component_data_view: 'typing.Optional[ComponentDataView]' = None

        if self.__value_name.find(":") == -1:
            # not found where the data location. So the default for value names
            # will be buffered (i.e. historical). So appending it
            self.__value_name += ":historical"

        if self.__tolerance < 0.0:
            raise RuntimeError(f"Tolerance should be non-negative [ tolerance = {self.__tolerance} ].")

        if self.__minimum_itr <= 0:
            raise RuntimeError(f"The minimum number of iterations needs to be greater than 0 [ minimum_itr = {self.__minimum_itr} ].")

        if self.__patience_itr >= self.__minimum_itr:
            raise RuntimeError(f"The number of patience iterations needs to be greater than minimum number of iterations [ patience_itr = {self.__patience_itr}, minimum_itr = {self.__minimum_itr} ].")

        self.__best_value: 'typing.Optional[float]' = None
        self.__patience_step = 0

    def Initialize(self):
        component = GetComponentHavingDataByFullName(self.__component_name, self.__optimization_problem)
        self.__component_data_view = ComponentDataView(component, self.__optimization_problem)

    @time_decorator()
    def IsConverged(self) -> bool:
        step = self.__optimization_problem.GetStep()
        value = GetComponentValueByFullName(self.__component_data_view, self.__value_name)

        if not isinstance(value, float) and not isinstance(value, int):
            raise RuntimeError(f"The value represented by {self.__value_name} is not a floating or int value.")

        if self.__best_value is None:
            self.__init_value = value
            self.__best_value = value

        if step >= self.__minimum_itr:
            if self.__best_value > value + self.__tolerance * self.__init_value:
                self.__best_value = value
                self.__patience_step = 0
            else:
                self.__patience_step += 1

            # say it is converged if you cannot find a better solution within the number of
            # patience iterations. That means, there were no improvement for patience iterations.
            self.__conv = self.__patience_step >= self.__patience_itr
        else:
            self.__conv = False

        self.__component_data_view.GetBufferedData().SetValue(f"{self.__value_name.split(':')[0]}_pat_itr_{self.__patience_itr}_best_value", self.__best_value)
        self.__component_data_view.GetBufferedData().SetValue(f"{self.__value_name.split(':')[0]}_pat_itr_{self.__patience_itr}_patience_step", self.__patience_step)

        return self.__conv

    def Finalize(self):
        pass

    def GetInfo(self) -> 'list[tuple[str, typing.Union[int, float, str]]]':
        info = [
                    ('type'          , 'patience_conv_criterion'),
                    ('best_value'    , self.__best_value),
                    ('patience_step' , self.__patience_step),
                    ('patience_itr'  , self.__patience_itr),
                    ('minimum_itr'   , self.__minimum_itr),
                    ('tolerance'     , self.__tolerance),
                    ('component_name', self.__component_name),
                    ('value_name'    , self.__value_name),
                    ('status'        , str("converged" if self.__conv else "not converged"))
                ]
        return info
