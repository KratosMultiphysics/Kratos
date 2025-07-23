import typing
import KratosMultiphysics as Kratos
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem import OptimizationProblem
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem_utilities import OptimizationComponentFactory
from KratosMultiphysics.OptimizationApplication.convergence_criterions.convergence_criteria import ConvergenceCriteria

def Factory(model: Kratos.Model, parameters: Kratos.Parameters, optimization_problem: OptimizationProblem) -> ConvergenceCriteria:
    if not parameters.Has("settings"):
        raise RuntimeError(f"CombinedConvCriteria instantiation requires a \"settings\" in parameters [ parameters = {parameters}].")
    return CombinedConvCriteria(model, parameters["settings"], optimization_problem)

class CombinedConvCriteria(ConvergenceCriteria):
    @classmethod
    def GetDefaultParameters(cls):
        return Kratos.Parameters("""{
            "operator"                      : "",
            "list_of_convergence_criterions": [
                {
                    "type"    : "",
                    "module"  : "KratosMultiphysics.OptimizationApplication.convergence_criterions",
                    "settings": {}
                }
            ]
        }""")

    def __init__(self, model: Kratos.Model, parameters: Kratos.Parameters, optimization_problem: OptimizationProblem):
        default_params = self.GetDefaultParameters()

        parameters.ValidateAndAssignDefaults(default_params)

        self.__operator = parameters["operator"].GetString()
        if self.__operator not in ["and", "or"]:
            raise RuntimeError(f"Unsupported operator = \"{self.__operator}\". Followings are supported:\n\tand\n\tor")

        self.__list_of_convergence_criterions: 'list[ConvergenceCriteria]' = []
        for sub_convergence_params in parameters["list_of_convergence_criterions"].values():
            sub_convergence_params.ValidateAndAssignDefaults(default_params["list_of_convergence_criterions"].values()[0])
            self.__list_of_convergence_criterions.append(OptimizationComponentFactory(model, sub_convergence_params, optimization_problem))

    def Add(self, convergence_criteria: ConvergenceCriteria) -> None:
        self.__list_of_convergence_criterions.append(convergence_criteria)

    def Initialize(self):
        if len(self.__list_of_convergence_criterions) == 0:
            raise RuntimeError("An empty combined convergence criteria is not allowed.")

        for conv in self.__list_of_convergence_criterions:
            conv.Initialize()

    def IsConverged(self) -> bool:
        if self.__operator == "and":
            return all([conv.IsConverged()  for conv in self.__list_of_convergence_criterions])
        elif self.__operator == "or":
            return any([conv.IsConverged()  for conv in self.__list_of_convergence_criterions])

    def Finalize(self):
        for conv in self.__list_of_convergence_criterions:
            conv.Finalize()

    def GetInfo(self) -> 'list[tuple[str, typing.Union[int, float, str]]]':
        if len(self.__list_of_convergence_criterions) > 1:
            info = [("(", "")]

        info.extend(self.__list_of_convergence_criterions[0].GetInfo())

        for conv in self.__list_of_convergence_criterions[1:]:
            info.append((f"--- {self.__operator} ---", ""))
            info.extend(conv.GetInfo())

        if len(self.__list_of_convergence_criterions) > 1:
            info = [(")", "")]

        return info

