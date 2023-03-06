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
            raise RuntimeError(f"No model parts were provided for LinearStrainEnergyResponseFunction.")

    def Check(self) -> None:
        # all of the model parts needs to have the same root model part
        root_model_part = self.model_parts[0].GetRootModelPart()
        for model_part in self.model_parts:
            if root_model_part != model_part.GetRootModelPart():
                raise RuntimeError(f"Root model part mismatch. Evaluated model parts must have the same root model part. [ Required root model part = {root_model_part.FullName()}, current model part = {model_part.FullName()} ]")

    def CalculateValue(self) -> float:
        # execute the primal analysis
        self.primal_analysis_execution_policy_wrapper.Execute()

        # computes the strain energy
        value = 0.0
        for model_part in self.model_parts:
            value += KratosOA.ResponseUtils.LinearStrainEnergyResponseUtils.CalculateValue(model_part)
        return value

    def CalculateSensitivity(self, sensitivity_variable: SupportedSensitivityFieldVariableTypes, sensitivity_model_part: Kratos.ModelPart) -> None:
        KratosOA.ResponseUtils.LinearStrainEnergyResponseUtils.CalculateSensitivity(sensitivity_model_part, self.perturbation_size, sensitivity_variable)


