from math import log
import KratosMultiphysics as Kratos
import KratosMultiphysics.OptimizationApplication as KratosOA
from KratosMultiphysics.OptimizationApplication.responses.response_function import ResponseFunction
from KratosMultiphysics.OptimizationApplication.responses.response_function import SupportedSensitivityFieldVariableTypes
from KratosMultiphysics.OptimizationApplication.utilities.union_utilities import SupportedSensitivityFieldVariableTypes
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem import OptimizationProblem
from KratosMultiphysics.OptimizationApplication.utilities.response_utilities import EvaluateValue
from KratosMultiphysics.OptimizationApplication.utilities.response_utilities import EvaluateGradient

class LogResponseFunction(ResponseFunction):
    def __init__(self, response_function: ResponseFunction, optimization_problem: OptimizationProblem):
        super().__init__(f"log({response_function.GetName()})")
        self.response_function = response_function
        self.optimization_problem = optimization_problem

    def GetImplementedPhysicalKratosVariables(self) -> 'list[SupportedSensitivityFieldVariableTypes]':
        return self.response_function.GetImplementedPhysicalKratosVariables()

    def Initialize(self) -> None:
        self.response_function.Initialize()

    def Check(self) -> None:
        self.response_function.Check()

    def Finalize(self) -> None:
        self.response_function.Finalize()

    def GetInfluencingModelPart(self) -> Kratos.ModelPart:
        return self.response_function.GetInfluencingModelPart()

    def CalculateValue(self) -> float:
        return log(EvaluateValue(self.response_function, self.optimization_problem))

    def CalculateGradient(self, physical_variable_collective_expressions: 'dict[SupportedSensitivityFieldVariableTypes, KratosOA.CollectiveExpression]') -> None:
        v = EvaluateValue(self.response_function, self.optimization_problem)
        resp_physical_variable_collective_expressions = EvaluateGradient(self.response_function, physical_variable_collective_expressions, self.optimization_problem)

        for variable, collective_expression in physical_variable_collective_expressions.items():
            for result, g in zip(collective_expression.GetContainerExpressions(), resp_physical_variable_collective_expressions[variable].GetContainerExpressions()):
                result.SetExpression((g / v).GetExpression())

    def GetChildResponses(self) -> 'list[ResponseFunction]':
        return [self.response_function]

    def __str__(self) -> str:
        return f"Response [type = {self.__class__.__name__}, name = {self.GetName()}]"