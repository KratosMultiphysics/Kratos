import KratosMultiphysics as Kratos
import KratosMultiphysics.OptimizationApplication as KratosOA
from KratosMultiphysics.OptimizationApplication.optimization_info import OptimizationInfo
from KratosMultiphysics.OptimizationApplication.responses.response_function import ResponseFunction
from KratosMultiphysics.OptimizationApplication.utilities.container_data import ContainerData


class DisplacementResponseFunction(ResponseFunction):
    def __init__(self, model: Kratos.Model, parameters: Kratos.Parameters, _: OptimizationInfo):
        default_settings = Kratos.Parameters("""{
            "evaluated_model_part_name": "PLEASE_PROVIDE_A_MODEL_PART_NAME"
        }""")
        parameters.ValidateAndAssignDefaults(default_settings)
        self.model_part = model[parameters["evaluated_model_part_name"].GetString()]

    def Check(self):
        data_communicator = self.model_part.GetCommunicator().GetDataCommunicator()
        if not KratosOA.OptimizationUtils.IsVariableExistsInAllContainerProperties(self.model_part.Elements, Kratos.YOUNG_MODULUS, data_communicator):
            raise RuntimeError(f"Some elements' properties in {self.model_part.FullName()} does not have YOUNG_MODULUS variable.")

        if not KratosOA.OptimizationUtils.AreAllEntitiesOfSameGeometryType(self.model_part.Elements, data_communicator):
            raise RuntimeError(
                f"{self.model_part.FullName()} has elements with different geometry types. Please break down this response to SumResponseFunction where each sub response function only has elements with one geometry type.")

    def CalculateValue(self) -> float:

        # execute the primal analysis
        self.primal_analysis_execution_policy_wrapper.Execute()

        print(self.model_part)
        # computes the strain energy
        return KratosOA.LinearStrainEnergyResponseUtilities.CalculateStrainEnergy(self.model_part)

        # print(f"{self}")

        # return KratosOA.MassResponseUtilities.CalculateMass(self.model_part)

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
