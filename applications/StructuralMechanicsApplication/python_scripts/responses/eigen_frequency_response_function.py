import KratosMultiphysics as Kratos
import KratosMultiphysics.StructuralMechanicsApplication as KratosStruct
from KratosMultiphysics.OptimizationApplication.responses.response_function import ResponseFunction
from KratosMultiphysics.OptimizationApplication.optimization_info import OptimizationInfo
from KratosMultiphysics.OptimizationApplication.utilities.execution_policy_decorator import ExecutionPolicyDecorator

class EigenFrequencyResponseFunction(ResponseFunction):
    def __init__(self, model: Kratos.Model, parameters: Kratos.Parameters, optimization_info: OptimizationInfo):
        super().__init__()

        default_parameters = Kratos.Parameters("""{
            "model_part_name"             : "PLEASE_PROVIDE_RESPONSE_MODEL_PART_NAME",
            "eigen_frequency_ids"         : [],
            "eigen_frequency_weights"     : [],
            "perturbation_size"           : 1e-8,
            "primal_execution_policy_name": "PLEASE_PROVIDE_PRIMAL_ANALYSIS_NAME"
        }""")

        parameters.ValidateAndAssignDefaults(default_parameters)

        self.model_part = model[parameters["model_part_name"].GetString()]
        self.primal_execution_policy: ExecutionPolicyDecorator = optimization_info.GetOptimizationProcess(ExecutionPolicyDecorator, parameters["primal_execution_policy_name"].GetString())

        perturbation_size = parameters["perturbation_size"].GetDouble()
        eigen_frequency_ids = [int(value) for value in parameters["eigen_frequency_ids"].GetVector()]
        eigen_frequency_weights = [float(value) for value in parameters["eigen_frequency_weights"].GetVector()]
        self.eigen_frequency_response_utils = KratosStruct.EigenfrequencyResponseFunctionUtility(self.model_part, perturbation_size, eigen_frequency_ids, eigen_frequency_weights)

    def Check(self):
        pass

    def CalculateValue(self) -> float:
        return self.eigen_frequency_response_utils.CalculateValue()

    def CalculateSensitivity(self, sensitivity_variable: any, sensitivity_model_part: Kratos.ModelPart):
        if sensitivity_model_part != self.model_part:
            raise RuntimeError(f"Response model part {self.model_part.FullName()} and sensitivity model part {sensitivity_model_part.FullName()} mismatch.")

        if sensitivity_variable == Kratos.SHAPE_SENSITIVITY:
            self.eigen_frequency_response_utils.CalculateGradient()
        elif sensitivity_variable.Name.endswith("_SENSITIVITY"):
            primal_variable = Kratos.KratosGlobals.GetVariable(sensitivity_variable.Name()[:-12])
            self.eigen_frequency_response_utils.CalculateEigenFrequencyMaterialVariableSensitivity(primal_variable, sensitivity_variable)
        else:
            raise RuntimeError(f"Unsupported sensitivity variable {sensitivity_variable.Name} requested for {sensitivity_model_part.FullName()} for {self.__class__.__name__} sensitivity computation.")

    def GetModelPart(self) -> Kratos.ModelPart:
        return self.model_part