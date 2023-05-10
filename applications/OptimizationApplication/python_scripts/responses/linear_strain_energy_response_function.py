import KratosMultiphysics as Kratos
import KratosMultiphysics.OptimizationApplication as KratosOA
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem import OptimizationProblem
from KratosMultiphysics.OptimizationApplication.execution_policies.execution_policy_decorator import ExecutionPolicyDecorator
from KratosMultiphysics.OptimizationApplication.responses.response_function import ResponseFunction
from KratosMultiphysics.OptimizationApplication.responses.response_function import SupportedSensitivityFieldVariableTypes


class LinearStrainEnergyResponseFunction(ResponseFunction):
    @classmethod
    def GetSensitivityFieldVariables(cls) -> 'list[SupportedSensitivityFieldVariableTypes]':
        return [Kratos.SHAPE_SENSITIVITY, KratosOA.YOUNG_MODULUS_SENSITIVITY, KratosOA.THICKNESS_SENSITIVITY, KratosOA.POISSON_RATIO_SENSITIVITY]

    def __init__(self, model: Kratos.Model, parameters: Kratos.Parameters, optimization_info: OptimizationProblem):
        default_settings = Kratos.Parameters("""{
            "evaluated_model_part_names": [
                "PLEASE_PROVIDE_A_MODEL_PART_NAME"
            ],
            "primal_analysis_name"      : "",
            "perturbation_size"         : 1e-8
        }""")
        parameters.ValidateAndAssignDefaults(default_settings)

        self.perturbation_size = parameters["perturbation_size"].GetDouble()
        self.primal_analysis_execution_policy_decorator: ExecutionPolicyDecorator = optimization_info.GetExecutionPolicy(parameters["primal_analysis_name"].GetString())

        self.model_parts: 'list[Kratos.ModelPart]' = []
        for model_part_name in parameters["evaluated_model_part_names"].GetStringArray():
            self.model_parts.append(model[model_part_name])

        if len(self.model_parts) == 0:
            raise RuntimeError("No model parts were provided for LinearStrainEnergyResponseFunction.")

    def Check(self) -> None:
        pass

    def CalculateValue(self) -> float:
        return KratosOA.ResponseUtils.LinearStrainEnergyResponseUtils.CalculateValue(self.model_parts)

    def CalculateSensitivity(self, sensitivity_model_part_variable_info: 'dict[SupportedSensitivityFieldVariableTypes, list[Kratos.ModelPart]]') -> None:
        # get the evaluated model part
        analysis_model_part = self.primal_analysis_execution_policy_decorator.GetExecutionPolicy().GetAnalysisModelPart()
        KratosOA.ResponseUtils.LinearStrainEnergyResponseUtils.CalculateSensitivity(analysis_model_part, self.model_parts, sensitivity_model_part_variable_info, self.perturbation_size)


