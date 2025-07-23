import typing
import KratosMultiphysics as Kratos
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem import OptimizationProblem
from KratosMultiphysics.OptimizationApplication.utilities.component_data_view import ComponentDataView
from KratosMultiphysics.OptimizationApplication.convergence_criteria.convergence_criterion import ConvergenceCriterion
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem_utilities import GetComponentHavingDataByFullName
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem_utilities import GetComponentValueByFullName
from KratosMultiphysics.OptimizationApplication.utilities.logger_utilities import time_decorator

def Factory(_: Kratos.Model, parameters: Kratos.Parameters, optimization_problem: OptimizationProblem) -> ConvergenceCriterion:
    if not parameters.Has("settings"):
        raise RuntimeError(f"ConstraintConvCriterion instantiation requires a \"settings\" in parameters [ parameters = {parameters}].")
    return ConstraintConvCriterion(parameters["settings"], optimization_problem)

class ConstraintConvCriterion(ConvergenceCriterion):
    """
    ConstraintConvCriterion is a convergence criterion for optimization problems, specifically designed to check the convergence of constraints.

    This class checks whether a specified value associated with a component in the optimization problem has reached a given tolerance.
    """
    @classmethod
    def GetDefaultParameters(cls):
        return Kratos.Parameters("""{
            "component_name": "",
            "value_name"    : "std_value",
            "tolerance"     : 0.0
        }""")

    def __init__(self, parameters: Kratos.Parameters, optimization_problem: OptimizationProblem):
        parameters.ValidateAndAssignDefaults(self.GetDefaultParameters())


        self.__component_name = parameters["component_name"].GetString()
        self.__value_name = parameters["value_name"].GetString()
        self.__tolerance = parameters["tolerance"].GetDouble()
        self.__optimization_problem = optimization_problem

        if self.__value_name.find(":") == -1:
            # not found where the data location. So the default for value names
            # will be buffered (i.e. historical). So appending it
            self.__value_name += ":historical"

        if self.__tolerance < 0.0:
            raise RuntimeError(f"Tolerance should be non-negative [ tolerance = {self.__tolerance} ].")

    def Initialize(self):
        component = GetComponentHavingDataByFullName(self.__component_name, self.__optimization_problem)
        self.__component_data_view = ComponentDataView(component, self.__optimization_problem)

    @time_decorator()
    def IsConverged(self) -> bool:
        self.__value = GetComponentValueByFullName(self.__component_data_view, self.__value_name)
        self.__conv = self.__value <= self.__tolerance
        return self.__conv

    def Finalize(self):
        pass

    def GetInfo(self) -> 'list[tuple[str, typing.Union[int, float, str]]]':
        info = [
                    ("type"          , "constraint_conv_criterion"),
                    ("value"         , self.__value),
                    ("tolerance"     , self.__tolerance),
                    ('component_name', self.__component_name),
                    ('value_name'    , self.__value_name),
                    ("status", str("converged" if self.__conv else "not converged"))
               ]
        return info