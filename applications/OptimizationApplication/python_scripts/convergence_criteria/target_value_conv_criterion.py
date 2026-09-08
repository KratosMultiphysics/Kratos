import typing
import KratosMultiphysics as Kratos
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem import OptimizationProblem
from KratosMultiphysics.OptimizationApplication.utilities.component_data_view import ComponentDataView
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem_utilities import GetComponentHavingDataByFullName
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem_utilities import GetComponentValueByFullName
from KratosMultiphysics.OptimizationApplication.convergence_criteria.convergence_criterion import ConvergenceCriterion
from KratosMultiphysics.OptimizationApplication.utilities.logger_utilities import time_decorator

def Factory(_: Kratos.Model, parameters: Kratos.Parameters, optimization_problem: OptimizationProblem) -> ConvergenceCriterion:
    if not parameters.Has("settings"):
        raise RuntimeError(f"TargetValueConvCriterion instantiation requires a \"settings\" in parameters [ parameters = {parameters}].")
    return TargetValueConvCriterion(parameters["settings"], optimization_problem)

class TargetValueConvCriterion(ConvergenceCriterion):
    """
    TargetValueConvCriterion is a convergence criterion for optimization problems that checks
    whether a specified value associated with a component has reached a target threshold.
    """
    @classmethod
    def GetDefaultParameters(cls):
        return Kratos.Parameters("""{
            "component_name": "algorithm",
            "value_name"    : "std_obj_value",
            "target_value"  : 1e-9
        }""")

    def __init__(self, parameters: Kratos.Parameters, optimization_problem: OptimizationProblem):
        parameters.ValidateAndAssignDefaults(self.GetDefaultParameters())

        self.__component_name = parameters["component_name"].GetString()
        self.__value_name = parameters["value_name"].GetString()
        self.__target_value = parameters["target_value"].GetDouble()
        self.__optimization_problem = optimization_problem
        self.__component_data_view: 'typing.Optional[ComponentDataView]' = None

        if self.__value_name.find(":") == -1:
            # not found where the data location. So the default for value names
            # will be buffered (i.e. historical). So appending it
            self.__value_name += ":historical"

    def Initialize(self):
        component = GetComponentHavingDataByFullName(self.__component_name, self.__optimization_problem)
        self.__component_data_view = ComponentDataView(component, self.__optimization_problem)

    @time_decorator()
    def IsConverged(self) -> bool:
        self.__current_value = GetComponentValueByFullName(self.__component_data_view, self.__value_name)
        self.__conv = self.__current_value <= self.__target_value
        return self.__conv

    def Finalize(self):
        pass

    def GetInfo(self) -> 'list[tuple[str, typing.Union[int, float, str]]]':
        info = [
                    ('type'          , 'target_value_conv_criteria'),
                    ('value'         , self.__current_value),
                    ('target_value'  , self.__target_value),
                    ('component_name', self.__component_name),
                    ('value_name'    , self.__value_name),
                    ('status'        , str("converged" if self.__conv else "not converged"))
               ]
        return info