from typing import Optional
from math import log

import KratosMultiphysics as Kratos
from KratosMultiphysics.OptimizationApplication.responses.response_function import ResponseFunction
from KratosMultiphysics.OptimizationApplication.responses.response_function import SupportedSensitivityFieldVariableTypes
from KratosMultiphysics.OptimizationApplication.utilities.model_part_utilities import ModelPartOperation
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem import OptimizationProblem
from KratosMultiphysics.OptimizationApplication.utilities.component_data_view import ComponentDataView
from KratosMultiphysics.OptimizationApplication.utilities.response_utilities import EvaluateValue
from KratosMultiphysics.OptimizationApplication.utilities.response_utilities import EvaluateGradient
from KratosMultiphysics.OptimizationApplication.utilities.response_utilities import BinaryOperator

class BinaryOperatorResponseFunction(ResponseFunction):
    def __init__(self, response_function_1: ResponseFunction, response_function_2: ResponseFunction, binary_operator: BinaryOperator, optimization_problem: OptimizationProblem):
        if binary_operator == BinaryOperator.DIVIDE:
            # this is because, the optimization_problem data container uses "/" as a path separator.
            super().__init__(f"({response_function_1.GetName()}÷{response_function_2.GetName()})")
        else:
            super().__init__(f"({response_function_1.GetName()}{binary_operator.value}{response_function_2.GetName()})")

        self.optimization_problem = optimization_problem
        self.response_function_1 = response_function_1
        self.response_function_2 = response_function_2
        self.binary_operator = binary_operator
        self.model_part: Optional[Kratos.ModelPart] = None

    def GetImplementedPhysicalKratosVariables(self) -> 'list[SupportedSensitivityFieldVariableTypes]':
        vars_list = self.response_function_1.GetImplementedPhysicalKratosVariables()
        vars_list.extend(self.response_function_2.GetImplementedPhysicalKratosVariables())
        return vars_list

    def Initialize(self) -> None:
        self.response_function_1.Initialize()
        self.response_function_2.Initialize()

        if len(self.response_function_1.GetImplementedPhysicalKratosVariables()) != 0 and len(self.response_function_2.GetImplementedPhysicalKratosVariables()) != 0:
            self.model_part = ModelPartOperation(self.response_function_1.GetInfluencingModelPart().GetModel(), ModelPartOperation.OperationType.UNION, f"response_{self.GetName()}", [self.response_function_1.GetInfluencingModelPart().FullName(), self.response_function_2.GetInfluencingModelPart().FullName()], False).GetModelPart()
        elif len(self.response_function_1.GetImplementedPhysicalKratosVariables()) != 0:
            self.model_part = self.response_function_1.GetInfluencingModelPart()
        elif len(self.response_function_2.GetImplementedPhysicalKratosVariables()) != 0:
            self.model_part = self.response_function_2.GetInfluencingModelPart()

    def Check(self) -> None:
        self.response_function_1.Check()
        self.response_function_2.Check()

    def Finalize(self) -> None:
        self.response_function_1.Finalize()
        self.response_function_2.Finalize()

    def GetInfluencingModelPart(self) -> Kratos.ModelPart:
        return self.model_part

    def CalculateValue(self) -> float:
        v1 = EvaluateValue(self.response_function_1, self.optimization_problem)
        v2 = EvaluateValue(self.response_function_2, self.optimization_problem)

        # now do the binary arithmetics.
        if self.binary_operator == BinaryOperator.ADD:
            return v1 + v2
        elif self.binary_operator == BinaryOperator.SUBTRACT:
            return v1 - v2
        elif self.binary_operator == BinaryOperator.MULTIPLY:
            return v1 * v2
        elif self.binary_operator == BinaryOperator.DIVIDE:
            return v1 / v2
        elif self.binary_operator == BinaryOperator.POWER:
            return v1 ** v2

    def CalculateGradient(self, physical_variable_combined_tensor_adaptor: 'dict[SupportedSensitivityFieldVariableTypes, Kratos.TensorAdaptors.DoubleCombinedTensorAdaptor]') -> None:
        v1 = EvaluateValue(self.response_function_1, self.optimization_problem)
        v2 = EvaluateValue(self.response_function_2, self.optimization_problem)

        resp_1_gradients = {variable: Kratos.TensorAdaptors.DoubleCombinedTensorAdaptor(cta) for variable, cta in physical_variable_combined_tensor_adaptor.items()}
        EvaluateGradient(self.response_function_1, resp_1_gradients, self.optimization_problem)

        resp_2_gradients = {variable: Kratos.TensorAdaptors.DoubleCombinedTensorAdaptor(cta) for variable, cta in physical_variable_combined_tensor_adaptor.items()}
        EvaluateGradient(self.response_function_2, resp_2_gradients, self.optimization_problem)

        for variable, result_cta in physical_variable_combined_tensor_adaptor.items():
            g1_cta = resp_1_gradients[variable]
            g2_cta = resp_2_gradients[variable]
            if self.binary_operator == BinaryOperator.ADD:
                result_cta.data = g1_cta.data + g2_cta.data
            elif self.binary_operator == BinaryOperator.SUBTRACT:
                result_cta.data = g1_cta.data - g2_cta.data
            elif self.binary_operator == BinaryOperator.MULTIPLY:
                result_cta.data = g1_cta.data * v2 + g2_cta.data * v1
            elif self.binary_operator == BinaryOperator.DIVIDE:
                result_cta.data = g1_cta.data / v2 - g2_cta.data * (v1 / v2 ** 2)
            elif self.binary_operator == BinaryOperator.POWER:
                result_cta.data = (g1_cta.data * (v2 / v1) + g2_cta.data * log(v1)) * (v1 ** v2)
            Kratos.TensorAdaptors.DoubleCombinedTensorAdaptor(result_cta, perform_store_data_recursively=False, copy=False).StoreData()

    def GetChildResponses(self) -> 'list[ResponseFunction]':
        return [self.response_function_1, self.response_function_2]

    def __str__(self) -> str:
        if self.model_part is not None:
            return f"Response [type = {self.__class__.__name__}, name = {self.GetName()}, model part name = {self.model_part.FullName()}]"
        else:
            return f"Response [type = {self.__class__.__name__}, name = {self.GetName()}, model part name = n/a ]"