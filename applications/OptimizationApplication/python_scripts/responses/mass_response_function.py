import KratosMultiphysics as Kratos
import KratosMultiphysics.OptimizationApplication as KratosOA
from KratosMultiphysics.OptimizationApplication.optimization_info import OptimizationInfo
from KratosMultiphysics.OptimizationApplication.responses.response_function import ResponseFunction
from KratosMultiphysics.OptimizationApplication.utilities.container_data import ContainerData


class MassResponseFunction(ResponseFunction):
    def __init__(self, model: Kratos.Model, parameters: Kratos.Parameters, _: OptimizationInfo):
        default_settings = Kratos.Parameters("""{
            "evaluated_model_part_name": "PLEASE_PROVIDE_A_MODEL_PART_NAME"
        }""")
        parameters.ValidateAndAssignDefaults(default_settings)
        self.model_part = model[parameters["evaluated_model_part_name"].GetString()]

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

    def CalculateSensitivity(self, sensitivity_variable, sensitivity_container: ContainerData):
        if sensitivity_variable == Kratos.SHAPE_SENSITIVITY and sensitivity_container.GetContainerTpe() == ContainerData.ContainerEnum.NODES:
            KratosOA.MassResponseUtilities.CalculateMassShapeSensitivity(sensitivity_container.GetModelPart(), sensitivity_variable)
        elif sensitivity_variable == KratosOA.DENSITY_SENSITIVITY and sensitivity_container.GetContainerTpe() == ContainerData.ContainerEnum.ELEMENT_PROPERTIES:
            KratosOA.MassResponseUtilities.CalculateMassDensitySensitivity(sensitivity_container.GetModelPart(), sensitivity_variable)
        elif sensitivity_variable == KratosOA.THICKNESS_SENSITIVITY and sensitivity_container.GetContainerTpe() == ContainerData.ContainerEnum.ELEMENT_PROPERTIES:
            KratosOA.MassResponseUtilities.CalculateMassThicknessSensitivity(sensitivity_container.GetModelPart(), sensitivity_variable)
        elif sensitivity_variable == KratosOA.CROSS_AREA_SENSITIVITY and sensitivity_container.GetContainerTpe() == ContainerData.ContainerEnum.ELEMENT_PROPERTIES:
            KratosOA.MassResponseUtilities.CalculateMassCrossAreaSensitivity(sensitivity_container.GetModelPart(), sensitivity_variable)
        else:
            msg = f"Unsupported sensitivity w.r.t. {sensitivity_variable.Name()} requested for {str(sensitivity_container)}."
            msg += "Followings are supported options:"
            msg += "\n\tSHAPE_SENSITIVITY for NODES container"
            msg += "\n\tDENSITY_SENSITIVITY for ELEMENT_PROPERTIES container"
            msg += "\n\tTHICKNESS_SENSITIVITY for ELEMENT_PROPERTIES container"
            msg += "\n\tCROSS_AREA_SENSITIVITY for ELEMENT_PROPERTIES container"
            raise RuntimeError(msg)

        # read the computed sensitivities to the vector
        sensitivity_container.ReadDataFromContainerVariable(sensitivity_variable)
