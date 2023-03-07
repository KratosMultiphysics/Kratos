import KratosMultiphysics as Kratos
import KratosMultiphysics.OptimizationApplication as KratosOA
from KratosMultiphysics.OptimizationApplication.utilities.optimization_info import OptimizationInfo
from KratosMultiphysics.OptimizationApplication.execution_policies.execution_policy_decorator import ExecutionPolicyDecorator
from KratosMultiphysics.OptimizationApplication.responses.response_function import ResponseFunction
from KratosMultiphysics.OptimizationApplication.responses.response_function import SupportedSensitivityFieldVariableTypes


class LinearStrainEnergyResponseFunction(ResponseFunction):
    @classmethod
    def GetSensitivityFieldVariables(cls) -> 'list[SupportedSensitivityFieldVariableTypes]':
        return [Kratos.SHAPE_SENSITIVITY, KratosOA.YOUNG_MODULUS_SENSITIVITY, KratosOA.POISSON_RATIO_SENSITIVITY]

    def __init__(self, model: Kratos.Model, parameters: Kratos.Parameters, optimization_info: OptimizationInfo):
        default_settings = Kratos.Parameters("""{
            "evaluated_model_part_names": [
                "PLEASE_PROVIDE_A_MODEL_PART_NAME"
            ],
            "primal_analysis_name"      : "",
            "perturbation_size"         : 1e-8
        }""")
        parameters.ValidateAndAssignDefaults(default_settings)

        self.primal_analysis_execution_policy_wrapper: ExecutionPolicyDecorator = optimization_info.GetOptimizationProcess(ExecutionPolicyDecorator, parameters["primal_analysis_name"].GetString())
        self.perturbation_size = parameters["perturbation_size"].GetDouble()

        self.model_parts: 'list[Kratos.ModelPart]' = []
        for model_part_name in parameters["evaluated_model_part_names"].GetStringArray():
            self.model_parts.append(model[model_part_name])

        if len(self.model_parts) == 0:
            raise RuntimeError("No model parts were provided for LinearStrainEnergyResponseFunction.")

    def Check(self) -> None:
        pass

    def CalculateValue(self) -> float:
        # execute the primal analysis
        self.primal_analysis_execution_policy_wrapper.Execute()

        # computes the strain energy
        return KratosOA.ResponseUtils.LinearStrainEnergyResponseUtils.CalculateValue(self.model_parts)

    def CalculateSensitivity(self, sensitivity_model_part_variable_info: 'dict[Kratos.ModelPart, list[SupportedSensitivityFieldVariableTypes]]') -> None:
        KratosOA.ResponseUtils.LinearStrainEnergyResponseUtils.CalculateSensitivity(self.model_parts, sensitivity_model_part_variable_info, self.perturbation_size)


