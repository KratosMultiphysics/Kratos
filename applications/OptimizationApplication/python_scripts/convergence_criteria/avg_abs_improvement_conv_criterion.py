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
        raise RuntimeError(f"AvgAbsImprovementConvCriterion instantiation requires a \"settings\" in parameters [ parameters = {parameters}].")
    return AvgAbsImprovementConvCriterion(parameters["settings"], optimization_problem)

class AvgAbsImprovementConvCriterion(ConvergenceCriterion):
    """
    AvgAbsImprovementConvCriterion is a convergence criterion for optimization algorithms that checks
    whether the average absolute improvement of a tracked value over a specified number of iterations
    falls below a given tolerance.
    """
    @classmethod
    def GetDefaultParameters(cls):
        return Kratos.Parameters("""{
            "component_name": "algorithm",
            "value_name"    : "std_obj_value",
            "tracked_iter"  : 5,
            "tolerance"     : 1e-6
        }""")

    def __init__(self, parameters: Kratos.Parameters, optimization_problem: OptimizationProblem):
        parameters.ValidateAndAssignDefaults(self.GetDefaultParameters())

        self.__component_name = parameters["component_name"].GetString()
        self.__value_name = parameters["value_name"].GetString()
        self.__tolerance = parameters["tolerance"].GetDouble()
        self.__tracked_iter = parameters["tracked_iter"].GetInt()
        self.__optimization_problem = optimization_problem
        self.__component_data_view: 'typing.Optional[ComponentDataView]' = None

        if self.__value_name.find(":") == -1:
            # not found where the data location. So the default for value names
            # will be buffered (i.e. historical). So appending it
            self.__value_name += ":historical"

        if self.__tolerance < 0.0:
            raise RuntimeError(f"Tolerance should be non-negative [ tolerance = {self.__tolerance} ].")

        if self.__tracked_iter <= 2:
            raise RuntimeError(f"The tracked_iter value must be greater than 2 [ tracked_iter = {self.__tracked_iter} ].")

        self.__abs_value_change = numpy.zeros(self.__tracked_iter - 1)

    def Initialize(self):
        component = GetComponentHavingDataByFullName(self.__component_name, self.__optimization_problem)
        self.__component_data_view = ComponentDataView(component, self.__optimization_problem)

    @time_decorator()
    def IsConverged(self) -> bool:
        step = self.__optimization_problem.GetStep()
        value = GetComponentValueByFullName(self.__component_data_view, self.__value_name)

        if not isinstance(value, float) and not isinstance(value, int):
            raise RuntimeError(f"The value represented by {self.__value_name} is not a floating or int value.")

        if step == 0:
            self.__init_value = value
            self.__value = 0.0
            self.__conv = False
        elif step > 0:
            prev_value = GetComponentValueByFullName(self.__component_data_view, self.__value_name, 1)
            self.__abs_value_change[step % (self.__tracked_iter - 1)] = ((value - prev_value) / self.__init_value)

            if step >= self.__tracked_iter - 1:
                self.__value = numpy.abs(numpy.sum(self.__abs_value_change)) / (self.__tracked_iter - 1)
                self.__conv = self.__value <= self.__tolerance
            else:
                self.__value = numpy.abs(numpy.sum(self.__abs_value_change)) / step
                self.__conv = False

        self.__component_data_view.GetBufferedData().SetValue(f"{self.__value_name.split(':')[0]}_avg_abs_{self.__tracked_iter}_itr", self.__value)

        return self.__conv

    def Finalize(self):
        pass

    def GetInfo(self) -> 'list[tuple[str, typing.Union[int, float, str]]]':
        info = [
                    ('type'          , 'avg_abs_improvement_conv_criterion'),
                    ('avg_value'     , self.__value),
                    ('tracked_iter'  , self.__tracked_iter),
                    ('tolerance'     , self.__tolerance),
                    ('component_name', self.__component_name),
                    ('value_name'    , self.__value_name),
                    ('status'        , str("converged" if self.__conv else "not converged"))
                ]
        return info
