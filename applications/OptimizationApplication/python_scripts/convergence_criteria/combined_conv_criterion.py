import typing
import KratosMultiphysics as Kratos
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem import OptimizationProblem
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem_utilities import OptimizationComponentFactory
from KratosMultiphysics.OptimizationApplication.convergence_criteria.convergence_criterion import ConvergenceCriterion

def Factory(model: Kratos.Model, parameters: Kratos.Parameters, optimization_problem: OptimizationProblem) -> ConvergenceCriterion:
    if not parameters.Has("settings"):
        raise RuntimeError(f"CombinedConvCriterion instantiation requires a \"settings\" in parameters [ parameters = {parameters}].")
    return CombinedConvCriterion(model, parameters["settings"], optimization_problem)

class CombinedConvCriterion(ConvergenceCriterion):
    """
    CombinedConvCriterion is a convergence criterion that combines multiple sub-criteria using a logical operator ("and" or "or").
    It evaluates convergence by applying the specified operator to the results of its sub-criteria.

    Args:
        model (Kratos.Model): The Kratos model instance.
        parameters (Kratos.Parameters): Parameters specifying the operator and the list of sub-criteria.
        optimization_problem (OptimizationProblem): The optimization problem instance.

    Raises:
        RuntimeError: If an unsupported operator is specified or if the list of sub-criteria is empty.
    """
    @classmethod
    def GetDefaultParameters(cls):
        return Kratos.Parameters("""{
            "operator"                    : "",
            "list_of_convergence_criteria": []
        }""")

    def __init__(self, model: Kratos.Model, parameters: Kratos.Parameters, optimization_problem: OptimizationProblem):
        default_params = self.GetDefaultParameters()

        parameters.ValidateAndAssignDefaults(default_params)

        self.__operator = parameters["operator"].GetString()
        if self.__operator not in ["and", "or"]:
            raise RuntimeError(f"Unsupported operator = \"{self.__operator}\". Followings are supported:\n\tand\n\tor")

        self.__list_of_convergence_criteria: 'list[ConvergenceCriterion]' = []
        default_sub_conv_settings = Kratos.Parameters("""{
            "type"    : "",
            "module"  : "KratosMultiphysics.OptimizationApplication.convergence_criteria",
            "settings": {}
        }""")
        for sub_convergence_params in parameters["list_of_convergence_criteria"].values():
            sub_convergence_params.ValidateAndAssignDefaults(default_sub_conv_settings)
            self.__list_of_convergence_criteria.append(OptimizationComponentFactory(model, sub_convergence_params, optimization_problem))

    def Add(self, convergence_criterion: ConvergenceCriterion) -> None:
        self.__list_of_convergence_criteria.append(convergence_criterion)

    def Initialize(self):
        if len(self.__list_of_convergence_criteria) == 0:
            raise RuntimeError("An empty combined convergence criteria is not allowed.")

        for conv in self.__list_of_convergence_criteria:
            conv.Initialize()

    def IsConverged(self) -> bool:
        if self.__operator == "and":
            self.__conv = all([conv.IsConverged()  for conv in self.__list_of_convergence_criteria])
        elif self.__operator == "or":
            self.__conv = any([conv.IsConverged()  for conv in self.__list_of_convergence_criteria])
        return self.__conv

    def Finalize(self):
        for conv in self.__list_of_convergence_criteria:
            conv.Finalize()

    def GetInfo(self) -> 'list[tuple[str, typing.Union[int, float, str]]]':
        info = []
        if len(self.__list_of_convergence_criteria) > 1:
            info.append(("(", ""))

        info.extend(self.__list_of_convergence_criteria[0].GetInfo())

        for conv in self.__list_of_convergence_criteria[1:]:
            info.append((f"--- {self.__operator} ---", ""))
            info.extend(conv.GetInfo())

        if len(self.__list_of_convergence_criteria) > 1:
            info.append((")  -->", str("converged" if self.__conv else "not converged")))

        return info

