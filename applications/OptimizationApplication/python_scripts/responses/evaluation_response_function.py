import KratosMultiphysics as Kratos
import KratosMultiphysics.OptimizationApplication as KratosOA
from KratosMultiphysics.OptimizationApplication.responses.response_function import ResponseFunction
from KratosMultiphysics.OptimizationApplication.responses.response_function import SupportedSensitivityFieldVariableTypes
from KratosMultiphysics.OptimizationApplication.utilities.union_utilities import SupportedSensitivityFieldVariableTypes
from KratosMultiphysics.OptimizationApplication.utilities.component_data_view import ComponentDataView
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem import OptimizationProblem
from KratosMultiphysics.OptimizationApplication.utilities.buffered_dict import BufferedDict

class EvaluationResponseFunction(ResponseFunction):
    def __init__(self, response_function: ResponseFunction, optimization_problem: OptimizationProblem):
        super().__init__(f"Eval_{response_function.GetName()}")
        self.response_function = response_function
        self.optimization_problem = optimization_problem

    def GetImplementedPhysicalKratosVariables(self) -> 'list[SupportedSensitivityFieldVariableTypes]':
        return self.response_function.GetImplementedPhysicalKratosVariables()

    def Initialize(self) -> None:
        self.response_function.Initialize()
        ComponentDataView("evaluated_responses", self.optimization_problem).SetDataBuffer(1)

    def Check(self) -> None:
        self.response_function.Check()

    def Finalize(self) -> None:
        self.response_function.Finalize()

    def GetInfluencingModelPart(self) -> Kratos.ModelPart:
        return self.response_function.GetInfluencingModelPart()

    def CalculateValue(self) -> float:
        # this response is used at the top most level of the evaluated chained responses.
        # so this creates a new data container under the optimization problem to avoid
        # having to compute the same response value twice.

        buffered_data = ComponentDataView("evaluated_responses", self.optimization_problem).GetBufferedData()

        # reset data of the evaluation
        self.__ResetEvaluationData(self, buffered_data, "values")

        # now calculate
        return self.response_function.CalculateValue()

    def CalculateGradient(self, physical_variable_collective_expressions: 'dict[SupportedSensitivityFieldVariableTypes, KratosOA.CollectiveExpression]') -> None:
        # this response is used at the top most level of the evaluated chained responses.
        # so this creates a new data container under the optimization problem to avoid
        # having to compute the same response gradient twice.

        unbuffered_data = ComponentDataView("evaluated_responses", self.optimization_problem).GetUnBufferedData()

        # reset data of the evaluation
        self.__ResetEvaluationData(self, unbuffered_data, "gradients")

        return self.response_function.CalculateGradient(physical_variable_collective_expressions)

    def GetChildResponses(self) -> 'list[ResponseFunction]':
        return [self.response_function]

    def __str__(self) -> str:
        return f"Response [type = {self.__class__.__name__}, name = {self.GetName()}]"

    @staticmethod
    def __ResetEvaluationData(response_function: ResponseFunction, buffered_data: BufferedDict, prefix: str) -> None:
        if buffered_data.HasValue(f"{prefix}/{response_function.GetName()}"):
            del buffered_data[f"{prefix}/{response_function.GetName()}"]
        for child_response in response_function.GetChildResponses():
            # now reset the children
            EvaluationResponseFunction.__ResetEvaluationData(child_response, buffered_data, prefix)
