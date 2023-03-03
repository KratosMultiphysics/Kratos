import KratosMultiphysics as Kratos
import KratosMultiphysics.OptimizationApplication as KratosOA
from KratosMultiphysics.OptimizationApplication.utilities.optimization_info import OptimizationInfo
from KratosMultiphysics.OptimizationApplication.execution_policies.execution_policy_decorator import ExecutionPolicyDecorator
from KratosMultiphysics.OptimizationApplication.responses.response_function import ResponseFunction


class LinearStrainEnergyResponseFunction(ResponseFunction):
    def __init__(self, model: Kratos.Model, parameters: Kratos.Parameters, optimization_info: OptimizationInfo):
        default_settings = Kratos.Parameters("""{
            "evaluated_model_part_name": "PLEASE_PROVIDE_A_MODEL_PART_NAME",
            "primal_analysis_name"     : "",
            "perturbation_size"        : 1e-8
        }""")
        parameters.ValidateAndAssignDefaults(default_settings)
        self.model_part = model[parameters["evaluated_model_part_name"].GetString()]
        self.primal_analysis_execution_policy_wrapper: ExecutionPolicyDecorator = optimization_info.GetOptimizationProcess(ExecutionPolicyDecorator, parameters["primal_analysis_name"].GetString())
        self.perturbation_size = parameters["perturbation_size"].GetDouble()

    def GetModelPart(self) -> Kratos.ModelPart:
        return self.model_part

    def Check(self) -> None:
        pass

    def CalculateValue(self) -> float:
        # execute the primal analysis
        self.primal_analysis_execution_policy_wrapper.Execute()

        # computes the strain energy
        return KratosOA.ResponseUtils.LinearStrainEnergyResponseUtils.CalculateValue(self.model_part)

    def CalculateSensitivity(self, sensitivity_variable: any, sensitivity_model_part: Kratos.ModelPart) -> None:
        KratosOA.ResponseUtils.LinearStrainEnergyResponseUtils.CalculateSensitivity(sensitivity_model_part, self.perturbation_size, sensitivity_variable)


