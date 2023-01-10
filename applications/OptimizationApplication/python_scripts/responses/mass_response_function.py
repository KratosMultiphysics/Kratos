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
from KratosMultiphysics.OptimizationApplication.responses.response_function import ResponseFunction
from KratosMultiphysics.OptimizationApplication.responses.response_function import ContainerEnum
from KratosMultiphysics.OptimizationApplication.responses.response_function import GetSensitivityContainer


class MassResponseFunction(ResponseFunction):
    def __init__(self, model: Kratos.Model, parameters: Kratos.Parameters):
        super().__init__(model, parameters)

        default_settings = Kratos.Parameters("""{
            "name"           : "",
            "model_part_name": "PLEASE_PROVIDE_A_MODEL_PART_NAME"
        }""")
        parameters.ValidateAndAssignDefaults(default_settings)
        self.model_part = self.model[parameters["model_part_name"].GetString()]
        self.name = parameters["name"].GetString()

    def Check(self):
        data_communicator = self.model_part.GetCommunicator().GetDataCommunicator()
        if not KratosOA.OptimizationVariableUtils.IsVariableExistsInAllContainerProperties(self.model_part.Elements, Kratos.DENSITY, data_communicator):
            raise RuntimeError(f"Some elements' properties in {self.model_part.FullName()} does not have DENSITY variable.")

        if KratosOA.OptimizationVariableUtils.IsVariableExistsInAtLeastOneContainerProperties(self.model_part.Elements, Kratos.THICKNESS, data_communicator) and \
           KratosOA.OptimizationVariableUtils.IsVariableExistsInAtLeastOneContainerProperties(self.model_part.Elements, Kratos.CROSS_AREA, data_communicator):
           raise RuntimeError(f"{self.model_part.FullName()} has elements consisting THICKNESS and CROSS_AREA. Please break down this response to SumResponseFunction where each sub response function only has elements with either THICKNESS or CROSS_AREA.")

        if not KratosOA.OptimizationVariableUtils.AreAllEntitiesOfSameGeometryType(self.model_part.Elements, data_communicator):
            raise RuntimeError(f"{self.model_part.FullName()} has elements with different geometry types. Please break down this response to SumResponseFunction where each sub response function only has elements with one geometry type.")

    def GetResponseFunctionName(self):
        return self.name

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
            Kratos.VariableUtils().SetNonHistoricalVariableToZero(sensitivity_variable, GetSensitivityContainer(sensitivity_container_type))


