import KratosMultiphysics as Kratos
import KratosMultiphysics.OptimizationApplication as KratosOA
from KratosMultiphysics.OptimizationApplication.responses.response_function import ResponseFunction

class MassResponseFunction(ResponseFunction):
    def __init__(self, model: Kratos.Model, parameters: Kratos.Parameters, _):
        default_settings = Kratos.Parameters("""{
            "evaluated_model_part_name": "PLEASE_PROVIDE_A_MODEL_PART_NAME"
        }""")
        parameters.ValidateAndAssignDefaults(default_settings)
        self.model_part = model[parameters["evaluated_model_part_name"].GetString()]

    def GetModelPart(self) -> Kratos.ModelPart:
        return self.model_part

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
        return KratosOA.ResponseUtils.MassResponseUtils.CalculateMass(self.model_part)

    def CalculateSensitivity(self, sensitivity_variable: any, sensitivity_model_part: Kratos.ModelPart):
        if sensitivity_variable == Kratos.SHAPE_SENSITIVITY:
            KratosOA.ResponseUtils.MassResponseUtils.CalculateMassShapeSensitivity(sensitivity_model_part, sensitivity_variable)
        elif sensitivity_variable == KratosOA.DENSITY_SENSITIVITY:
            KratosOA.ResponseUtils.MassResponseUtils.CalculateMassDensitySensitivity(sensitivity_model_part, sensitivity_variable)
        elif sensitivity_variable == KratosOA.THICKNESS_SENSITIVITY:
            KratosOA.ResponseUtils.MassResponseUtils.CalculateMassThicknessSensitivity(sensitivity_model_part, sensitivity_variable)
        elif sensitivity_variable == KratosOA.CROSS_AREA_SENSITIVITY:
            KratosOA.ResponseUtils.MassResponseUtils.CalculateMassCrossAreaSensitivity(sensitivity_model_part, sensitivity_variable)
        elif sensitivity_variable in [KratosOA.YOUNG_MODULUS_SENSITIVITY]:
            zero_sensitivity = KratosOA.ElementPropertiesContainerVariableDataHolder(sensitivity_model_part)
            zero_sensitivity.SetDataForContainerVariableToZero(KratosOA.YOUNG_MODULUS_SENSITIVITY)
            zero_sensitivity.AssignDataToContainerVariable(KratosOA.YOUNG_MODULUS_SENSITIVITY)
        else:
            msg = f"Unsupported sensitivity w.r.t. {sensitivity_variable.Name()} requested for {sensitivity_model_part.FullName()}. "
            msg += "Followings are supported variables:"
            msg += "\n\tSHAPE_SENSITIVITY"
            msg += "\n\tDENSITY_SENSITIVITY"
            msg += "\n\tTHICKNESS_SENSITIVITY"
            msg += "\n\tCROSS_AREA_SENSITIVITY"
            msg += "\n\tYOUNG_MODULUS_SENSITIVITY"
            raise RuntimeError(msg)

