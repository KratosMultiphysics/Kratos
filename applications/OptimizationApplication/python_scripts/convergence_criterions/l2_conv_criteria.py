import typing
import numpy
import KratosMultiphysics as Kratos
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem import OptimizationProblem
from KratosMultiphysics.OptimizationApplication.utilities.component_data_view import ComponentDataView
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem_utilities import GetComponentHavingDataByFullName
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem_utilities import GetComponentValueByFullName
from KratosMultiphysics.OptimizationApplication.convergence_criterions.convergence_criteria import ConvergenceCriteria
from KratosMultiphysics.OptimizationApplication.utilities.logger_utilities import time_decorator

def Factory(_: Kratos.Model, parameters: Kratos.Parameters, optimization_problem: OptimizationProblem) -> ConvergenceCriteria:
    if not parameters.Has("settings"):
        raise RuntimeError(f"L2ConvCriteria instantiation requires a \"settings\" in parameters [ parameters = {parameters}].")
    return L2ConvCriteria(parameters["settings"], optimization_problem)

class L2ConvCriteria(ConvergenceCriteria):
    @classmethod
    def GetDefaultParameters(cls):
        return Kratos.Parameters("""{
            "component_name": "algorithm",
            "field_name"    : "search_direction",
            "tolerance"     : 1e-9
        }""")

    def __init__(self, parameters: Kratos.Parameters, optimization_problem: OptimizationProblem):
        parameters.ValidateAndAssignDefaults(self.GetDefaultParameters())

        self.__component_name = parameters["component_name"].GetString()
        self.__field_name = parameters["field_name"].GetString()
        self.__tolerance = parameters["tolerance"].GetDouble()
        self.__optimization_problem = optimization_problem
        self.__component_data_view: 'typing.Optional[ComponentDataView]' = None

        if self.__field_name.find(":") == -1:
            # not found where the data location. So the default for field names
            # will be unbuffered (i.e. non_historical). So appending it
            self.__field_name += ":non_historical"

        if self.__tolerance < 0.0:
            raise RuntimeError(f"Tolerance should be non-negative [ tolerance = {self.__tolerance} ].")

    def Initialize(self):
        component = GetComponentHavingDataByFullName(self.__component_name, self.__optimization_problem)
        self.__component_data_view = ComponentDataView(component, self.__optimization_problem)

    @time_decorator()
    def IsConverged(self) -> bool:
        field = GetComponentValueByFullName(self.__component_data_view, self.__field_name)

        if not hasattr(field, "Evaluate"):
            raise RuntimeError(f"The value represented by {self.__field_name} is not a field.")

        self.__norm = numpy.linalg.norm(field.Evaluate().flatten())
        self.__conv = self.__norm <= self.__tolerance
        self.__component_data_view.GetBufferedData().SetValue(self.__field_name.split(":")[0] + "_l2_norm", self.__norm)
        return self.__conv

    def Finalize(self):
        pass

    def GetInfo(self) -> 'list[tuple[str, typing.Union[int, float, str]]]':
        info = [
                    ('type'          , 'l2_conv_criteria'),
                    ('l2_norm'       , self.__norm),
                    ('tolerance'     , self.__tolerance),
                    ('component_name', self.__component_name),
                    ('field_name'    , self.__field_name),
                    ('status'        , str("converged" if self.__conv else "not converged"))
               ]
        return info
