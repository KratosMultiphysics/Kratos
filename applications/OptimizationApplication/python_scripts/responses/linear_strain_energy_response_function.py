# ==============================================================================
#  KratosOptimizationApplication
#
#  License:         BSD License
#                   license: OptimizationApplication/license.txt
#
#  Main authors:    Suneth Warnakulasuriya
#
# ==============================================================================

import KratosMultiphysics as Kratos
import KratosMultiphysics.OptimizationApplication as KratosOA
from KratosMultiphysics.OptimizationApplication.optimization_info import OptimizationInfo
from KratosMultiphysics.OptimizationApplication.execution_policies.execution_policy_wrapper import ExecutionPolicyWrapper
from KratosMultiphysics.OptimizationApplication.responses.response_function import ResponseFunction
from KratosMultiphysics.OptimizationApplication.utilities.helper_utils import ContainerEnum
from KratosMultiphysics.OptimizationApplication.utilities.helper_utils import GetSensitivityContainer


class LinearStrainEnergyResponseFunction(ResponseFunction):
    def __init__(self, model: Kratos.Model, parameters: Kratos.Parameters, optimization_info: OptimizationInfo):
        super().__init__(model, parameters, optimization_info)

        default_settings = Kratos.Parameters("""{
            "model_part_name"      : "PLEASE_PROVIDE_A_MODEL_PART_NAME",
            "primal_analysis_name" : "",
            "perturbation_size"    : 1e-8
        }""")
        parameters.ValidateAndAssignDefaults(default_settings)
        self.model_part = self.model[parameters["model_part_name"].GetString()]
        self.primal_analysis_execution_policy_wrapper: ExecutionPolicyWrapper = optimization_info.GetRoutine("ExecutionPolicyWrapper", parameters["primal_analysis_name"].GetString())
        self.perturbation_size = parameters["perturbation_size"].GetDouble()

    def Check(self):
        pass

    def CalculateValue(self) -> float:
        # execute the primal analysis
        self.primal_analysis_execution_policy_wrapper.Execute()

        # computes the strain energy
        return KratosOA.LinearStrainEnergyResponseUtilities.CalculateStrainEnergy(self.model_part)

    def CalculateSensitivity(self, sensitivity_variable, sensitivity_model_part: Kratos.ModelPart, sensitivity_container_type: ContainerEnum) -> None:
        # computes strain energy derivatives
        if sensitivity_variable == Kratos.SHAPE_SENSITIVITY and sensitivity_container_type == ContainerEnum.NODES:
            KratosOA.LinearStrainEnergyResponseUtilities.CalculateStrainEnergyShapeSensitivity(sensitivity_model_part, self.perturbation_size, sensitivity_variable)
        elif sensitivity_container_type == ContainerEnum.ELEMENT_PROPERTIES:
            variable_name: str = sensitivity_variable.Name()
            if variable_name.endswith("_SENSITIVITY"):
                primal_variable = Kratos.KratosGlobals.GetVariable(variable_name[:-12])
                KratosOA.LinearStrainEnergyResponseUtilities.CalculateStrainEnergyElementPropertiesSensitivity(sensitivity_model_part, self.perturbation_size, primal_variable, sensitivity_variable)
            else:
                raise RuntimeError(f"{variable_name} is not a supported sensitivity variable.")
        else:
            # since this response function covers all the possible sensitivity calculations for the method given in CalculateMass,
            # we can return the corresponding zero values for other sensitivity values
            Kratos.VariableUtils().SetNonHistoricalVariableToZero(sensitivity_variable, GetSensitivityContainer(sensitivity_container_type))


