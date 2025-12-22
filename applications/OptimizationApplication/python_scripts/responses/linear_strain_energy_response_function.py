from typing import Optional

import KratosMultiphysics as Kratos
import KratosMultiphysics.OptimizationApplication as KratosOA
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem import OptimizationProblem
from KratosMultiphysics.OptimizationApplication.execution_policies.execution_policy_decorator import ExecutionPolicyDecorator
from KratosMultiphysics.OptimizationApplication.responses.response_function import ResponseFunction
from KratosMultiphysics.OptimizationApplication.responses.response_function import SupportedSensitivityFieldVariableTypes
from KratosMultiphysics.OptimizationApplication.utilities.model_part_utilities import ModelPartOperation
from KratosMultiphysics.OptimizationApplication.utilities.model_part_utilities import ModelPartUtilities

def Factory(model: Kratos.Model, parameters: Kratos.Parameters, optimization_problem: OptimizationProblem) -> ResponseFunction:
    if not parameters.Has("name"):
        raise RuntimeError(f"LinearStrainEnergyResponseFunction instantiation requires a \"name\" in parameters [ parameters = {parameters}].")
    if not parameters.Has("settings"):
        raise RuntimeError(f"LinearStrainEnergyResponseFunction instantiation requires a \"settings\" in parameters [ parameters = {parameters}].")

    return LinearStrainEnergyResponseFunction(parameters["name"].GetString(), model, parameters["settings"], optimization_problem)

class LinearStrainEnergyResponseFunction(ResponseFunction):
    def __init__(self, name: str, model: Kratos.Model, parameters: Kratos.Parameters, optimization_problem: OptimizationProblem):
        super().__init__(name)

        default_settings = Kratos.Parameters("""{
            "primal_analysis_name"           : "",
            "perturbation_size"              : 1e-8,
            "evaluated_model_part_names"     : [
                "PLEASE_PROVIDE_A_MODEL_PART_NAME"
            ]
        }""")
        parameters.ValidateAndAssignDefaults(default_settings)

        self.perturbation_size = parameters["perturbation_size"].GetDouble()
        self.model = model
        self.primal_analysis_execution_policy_decorator: ExecutionPolicyDecorator = optimization_problem.GetExecutionPolicy(parameters["primal_analysis_name"].GetString())

        evaluated_model_part_names = parameters["evaluated_model_part_names"].GetStringArray()
        if len(evaluated_model_part_names) == 0:
            raise RuntimeError(f"No model parts were provided for LinearStrainEnergyResponseFunction. [ response name = \"{self.GetName()}\"]")

        self.model_part_operation = ModelPartOperation(self.model, ModelPartOperation.OperationType.UNION, f"response_{self.GetName()}", evaluated_model_part_names, False)
        self.model_part: Optional[Kratos.ModelPart] = None

    def GetImplementedPhysicalKratosVariables(self) -> 'list[SupportedSensitivityFieldVariableTypes]':
        return [KratosOA.SHAPE, Kratos.YOUNG_MODULUS, Kratos.THICKNESS, Kratos.POISSON_RATIO]

    def Initialize(self) -> None:
        self.model_part = self.model_part_operation.GetModelPart()

    def Check(self) -> None:
        pass

    def Finalize(self) -> None:
        pass

    def GetInfluencingModelPart(self) -> Kratos.ModelPart:
        return self.primal_analysis_execution_policy_decorator.GetAnalysisModelPart()

    def CalculateValue(self) -> float:
        self.primal_analysis_execution_policy_decorator.Execute()
        return KratosOA.ResponseUtils.LinearStrainEnergyResponseUtils.CalculateValue(self.model_part)

    def CalculateGradient(self, physical_variable_collective_expressions: 'dict[SupportedSensitivityFieldVariableTypes, KratosOA.CollectiveExpression]') -> None:
        # first merge all the model parts
        merged_model_part_map = ModelPartUtilities.GetMergedMap(physical_variable_collective_expressions, False)

        # now get the intersected model parts
        intersected_model_part_map = ModelPartUtilities.GetIntersectedMap(self.model_part, merged_model_part_map, True)

        # TODO: Remove this block once everything is updated to TensorAdaptors
        # =========================================================
        physical_variable_combined_tensor_adaptors: 'dict[SupportedSensitivityFieldVariableTypes, Kratos.TensorAdaptors.DoubleCombinedTensorAdaptor]' = {}
        for physical_variable, collective_expression in physical_variable_collective_expressions.items():
            ta_list = []
            for ce in collective_expression.GetContainerExpressions():
                if isinstance(physical_variable, Kratos.DoubleVariable):
                    shape = Kratos.DenseVectorUnsignedInt([len(ce.GetContainer())])
                elif isinstance(physical_variable, Kratos.Array1DVariable3):
                    shape = Kratos.DenseVectorUnsignedInt([len(ce.GetContainer()), 3])
                else:
                    raise RuntimeError("Unsupported type")

                ta_list.append(Kratos.TensorAdaptors.DoubleTensorAdaptor(ce.GetContainer(), Kratos.DoubleNDData(shape), copy=False))
            cta = Kratos.TensorAdaptors.DoubleCombinedTensorAdaptor(ta_list, False, False)
            cta.CollectData()
            physical_variable_combined_tensor_adaptors[physical_variable] = cta
        # =========================================================

        # calculate the gradients
        for physical_variable, merged_model_part in merged_model_part_map.items():
            KratosOA.ResponseUtils.LinearStrainEnergyResponseUtils.CalculateGradient(
                physical_variable,
                merged_model_part,
                intersected_model_part_map[physical_variable],
                physical_variable_combined_tensor_adaptors[physical_variable],
                self.perturbation_size)

        # TODO: Remove this block once everything is updated to TensorAdaptors
        # =========================================================
        for physical_variable, collective_expression in physical_variable_collective_expressions.items():
            for i, ce in enumerate(collective_expression.GetContainerExpressions()):
                Kratos.Expression.CArrayExpressionIO.Read(ce, physical_variable_combined_tensor_adaptors[physical_variable].GetTensorAdaptors()[i].data)
        # =========================================================

    def __str__(self) -> str:
        return f"Response [type = {self.__class__.__name__}, name = {self.GetName()}, model part name = {self.model_part.FullName()}]"