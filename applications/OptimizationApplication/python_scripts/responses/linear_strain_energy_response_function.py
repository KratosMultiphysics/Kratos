import KratosMultiphysics as Kratos
import KratosMultiphysics.OptimizationApplication as KratosOA
from KratosMultiphysics.OptimizationApplication.optimization_info import OptimizationInfo
from KratosMultiphysics.OptimizationApplication.execution_policies.execution_policy_wrapper import ExecutionPolicyWrapper
from KratosMultiphysics.OptimizationApplication.responses.response_function import ResponseFunction
from KratosMultiphysics.OptimizationApplication.utilities.container_data import ContainerData


class LinearStrainEnergyResponseFunction(ResponseFunction):
    def __init__(self, model: Kratos.Model, parameters: Kratos.Parameters, optimization_info: OptimizationInfo):
        default_settings = Kratos.Parameters("""{
            "evaluated_model_part_name": "PLEASE_PROVIDE_A_MODEL_PART_NAME",
            "primal_analysis_name"     : "",
            "perturbation_size"        : 1e-8
        }""")
        parameters.ValidateAndAssignDefaults(default_settings)
        self.model_part = model[parameters["evaluated_model_part_name"].GetString()]
        self.primal_analysis_execution_policy_wrapper: ExecutionPolicyWrapper = optimization_info.GetOptimizationRoutine(ExecutionPolicyWrapper, parameters["primal_analysis_name"].GetString())
        self.perturbation_size = parameters["perturbation_size"].GetDouble()

    def Check(self):
        pass

    def CalculateValue(self) -> float:
        # execute the primal analysis
        self.primal_analysis_execution_policy_wrapper.Execute()

        # computes the strain energy
        return KratosOA.LinearStrainEnergyResponseUtilities.CalculateStrainEnergy(self.model_part)

    def CalculateSensitivity(self, sensitivity_variable, sensitivity_container: ContainerData):
        self.primal_analysis_execution_policy_wrapper.Execute()

        # computes strain energy derivatives
        if sensitivity_variable == Kratos.SHAPE_SENSITIVITY and sensitivity_container.GetContainerTpe() == ContainerData.ContainerEnum.NODES:
            KratosOA.LinearStrainEnergyResponseUtilities.CalculateStrainEnergyShapeSensitivity(sensitivity_container.GetModelPart(), self.perturbation_size, sensitivity_variable)
        elif sensitivity_container.GetContainerTpe() == ContainerData.ContainerEnum.ELEMENT_PROPERTIES:
            variable_name: str = sensitivity_variable.Name()
            if variable_name.endswith("_SENSITIVITY"):
                primal_variable = Kratos.KratosGlobals.GetVariable(variable_name[:-12])
                KratosOA.LinearStrainEnergyResponseUtilities.CalculateStrainEnergyElementPropertiesSensitivity(sensitivity_container.GetModelPart(), self.perturbation_size, primal_variable, sensitivity_variable)
            else:
                raise RuntimeError(f"{variable_name} is not a supported sensitivity variable.")
        else:
            msg = f"Unsupported sensitivity w.r.t. {sensitivity_variable.Name()} requested for {str(sensitivity_container)}."
            msg += "Followings are supported options:"
            msg += "\n\tSHAPE_SENSITIVITY for NODES container"
            msg += "\n\tmaterial properties for ELEMENT_PROPERTIES container"
            raise RuntimeError(msg)

        sensitivity_container.ReadDataFromContainer(sensitivity_variable)


