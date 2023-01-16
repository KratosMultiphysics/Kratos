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
from KratosMultiphysics.OptimizationApplication.responses.response_function import ResponseFunction
from KratosMultiphysics.OptimizationApplication.utilities.helper_utils import ContainerEnum
from KratosMultiphysics.OptimizationApplication.utilities.helper_utils import GetSensitivityContainer


class MassResponseFunction(ResponseFunction):
    def __init__(self, model: Kratos.Model, parameters: Kratos.Parameters, optimization_info: OptimizationInfo):
        super().__init__(model, parameters, optimization_info)

        default_settings = Kratos.Parameters("""{
            "model_part_name": "PLEASE_PROVIDE_A_MODEL_PART_NAME"
        }""")
        parameters.ValidateAndAssignDefaults(default_settings)
        self.model_part = self.model[parameters["model_part_name"].GetString()]

    def Check(self):
        data_communicator = self.model_part.GetCommunicator().GetDataCommunicator()
        if not KratosOA.OptimizationUtils.IsVariableExistsInAllContainerProperties(self.model_part.Elements, Kratos.DENSITY, data_communicator):
            raise RuntimeError(f"Some elements' properties in {self.model_part.FullName()} does not have DENSITY variable.")

        if KratosOA.OptimizationUtils.IsVariableExistsInAtLeastOneContainerProperties(self.model_part.Elements, Kratos.THICKNESS, data_communicator) and \
           KratosOA.OptimizationUtils.IsVariableExistsInAtLeastOneContainerProperties(self.model_part.Elements, KratosOA.CROSS_AREA, data_communicator):
           raise RuntimeError(f"{self.model_part.FullName()} has elements consisting THICKNESS and CROSS_AREA. Please break down this response to SumResponseFunction where each sub response function only has elements with either THICKNESS or CROSS_AREA.")

        if not KratosOA.OptimizationUtils.AreAllEntitiesOfSameGeometryType(self.model_part.Elements, data_communicator):
            raise RuntimeError(f"{self.model_part.FullName()} has elements with different geometry types. Please break down this response to SumResponseFunction where each sub response function only has elements with one geometry type.")

    def CalculateValue(self) -> float:
        return KratosOA.MassResponseUtilities.CalculateMass(self.model_part)

    def CalculateSensitivity(self, sensitivity_variable, sensitivity_model_part: Kratos.ModelPart, sensitivity_container_type: ContainerEnum) -> None:
        if sensitivity_variable == Kratos.SHAPE_SENSITIVITY and sensitivity_container_type == ContainerEnum.NODES:
            KratosOA.MassResponseUtilities.CalculateMassShapeSensitivity(sensitivity_model_part, sensitivity_variable)
        elif sensitivity_variable == KratosOA.DENSITY_SENSITIVITY and sensitivity_container_type == ContainerEnum.ELEMENT_PROPERTIES:
            KratosOA.MassResponseUtilities.CalculateMassDensitySensitivity(sensitivity_model_part, sensitivity_variable)
        elif sensitivity_variable == KratosOA.THICKNESS_SENSITIVITY and sensitivity_container_type == ContainerEnum.ELEMENT_PROPERTIES:
            KratosOA.MassResponseUtilities.CalculateMassThicknessSensitivity(sensitivity_model_part, sensitivity_variable)
        elif sensitivity_variable == KratosOA.CROSS_AREA_SENSITIVITY and sensitivity_container_type == ContainerEnum.ELEMENT_PROPERTIES:
            KratosOA.MassResponseUtilities.CalculateMassCrossAreaSensitivity(sensitivity_model_part, sensitivity_variable)
        else:
            # since this response function covers all the possible sensitivity calculations for the method given in CalculateMass,
            # we can return the corresponding zero values for other sensitivity values
            Kratos.VariableUtils().SetNonHistoricalVariableToZero(sensitivity_variable, GetSensitivityContainer(sensitivity_model_part, sensitivity_container_type))


