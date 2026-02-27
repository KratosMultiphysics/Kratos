from math import log
import KratosMultiphysics as Kratos
from KratosMultiphysics.OptimizationApplication.responses.response_function import ResponseFunction
from KratosMultiphysics.OptimizationApplication.responses.response_function import SupportedSensitivityFieldVariableTypes
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

    def CalculateGradient(self, physical_variable_combined_tensor_adaptor: 'dict[SupportedSensitivityFieldVariableTypes, Kratos.TensorAdaptors.DoubleCombinedTensorAdaptor]') -> None:
        v = EvaluateValue(self.response_function, self.optimization_problem)

        resp_gradients = {var: Kratos.TensorAdaptors.DoubleCombinedTensorAdaptor(cta) for var, cta in physical_variable_combined_tensor_adaptor.items()}
        EvaluateGradient(self.response_function, resp_gradients, self.optimization_problem)

        for variable, cta in physical_variable_combined_tensor_adaptor.items():
            cta.data = resp_gradients[variable].data / v
            Kratos.TensorAdaptors.DoubleCombinedTensorAdaptor(cta, perform_collect_data_recursively=False, perform_store_data_recursively=False, copy=False).StoreData()

    def GetChildResponses(self) -> 'list[ResponseFunction]':
        return [self.response_function]

    def __str__(self) -> str:
        return f"Response [type = {self.__class__.__name__}, name = {self.GetName()}]"