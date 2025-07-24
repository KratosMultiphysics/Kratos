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
        raise RuntimeError(f"MagnitudeReductionConvCriterion instantiation requires a \"settings\" in parameters [ parameters = {parameters}].")
    return MagnitudeReductionConvCriterion(parameters["settings"], optimization_problem)

class MagnitudeReductionConvCriterion(ConvergenceCriterion):
    """
    Convergence criterion based on the reduction of a magnitude (e.g., objective value) during optimization.

    This criterion checks whether the current value of a specified component has been reduced below a target value,
    which is computed as a scaling factor of the initial value or a specified machine precision, whichever is greater.

    Raises:
        RuntimeError: If the target scaling factor is not greater than zero or machine precision is negative.
    """
    @classmethod
    def GetDefaultParameters(cls):
        return Kratos.Parameters("""{
            "component_name"       : "algorithm",
            "value_name"           : "std_obj_value",
            "target_scaling_factor": 1e-3,
            "machine_precision"    : 1e-16
        }""")

    def __init__(self, parameters: Kratos.Parameters, optimization_problem: OptimizationProblem):
        parameters.ValidateAndAssignDefaults(self.GetDefaultParameters())

        self.__component_name = parameters["component_name"].GetString()
        self.__value_name = parameters["value_name"].GetString()
        self.__target_scaling_factor = parameters["target_scaling_factor"].GetDouble()
        self.__machine_precision = parameters["machine_precision"].GetDouble()
        self.__optimization_problem = optimization_problem
        self.__component_data_view: 'typing.Optional[ComponentDataView]' = None

        if self.__value_name.find(":") == -1:
            # not found where the data location. So the default for value names
            # will be buffered (i.e. historical). So appending it
            self.__value_name += ":historical"

        if self.__target_scaling_factor <= 0.0:
            raise RuntimeError(f"The target scaling factor should be greater than zero [ target_scaling_factor = {self.__target_scaling_factor} ].")

        if self.__machine_precision < 0.0:
            raise RuntimeError(f"The machine precision should be non-negative [ machine_precision = {self.__machine_precision} ].")

    def Initialize(self):
        component = GetComponentHavingDataByFullName(self.__component_name, self.__optimization_problem)
        self.__component_data_view = ComponentDataView(component, self.__optimization_problem)

    @time_decorator()
    def IsConverged(self) -> bool:
        iter = self.__optimization_problem.GetStep()
        self.__current_value = GetComponentValueByFullName(self.__component_data_view, self.__value_name)

        if iter == 0:
            self.__target_value = max(self.__current_value * self.__target_scaling_factor, self.__machine_precision)
            self.__component_data_view.GetUnBufferedData().SetValue(f"{self.__value_name.split(':')[0]}_conv_target", self.__target_value)
            self.__conv = False
        else:
            self.__conv = self.__current_value <= self.__target_value

        return self.__conv

    def Finalize(self):
        pass

    def GetInfo(self) -> 'list[tuple[str, typing.Union[int, float, str]]]':
        info = [
                    ('type'          , 'magnitude_reduction_conv_criterion'),
                    ('value'         , self.__current_value),
                    ('target_value'  , self.__target_value),
                    ('component_name', self.__component_name),
                    ('value_name'    , self.__value_name),
                    ('status'        , str("converged" if self.__conv else "not converged"))
               ]
        return info