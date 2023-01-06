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
from KratosMultiphysics.OptimizationApplication.responses.response_function import ResponesFunction
from KratosMultiphysics.OptimizationApplication.responses.response_sensitivity import ResponseSentivity
from KratosMultiphysics.OptimizationApplication.responses.response_sensitivity import ContainerEnum

class MassResponseFunction(ResponesFunction):
    def __init__(self, model: Kratos.Model, parameters: Kratos.Parameters):
        super().__init__(model, parameters)

        default_settings = Kratos.Parameters("""{
            "name"                      : "",
            "model_part_name"           : "PLEASE_PROVIDE_A_MODEL_PART_NAME",
            "sensitivity_variable_names": [""]
        }""")
        parameters.ValidateAndAssignDefaults(default_settings)
        self.model_part = self.model[parameters["model_part_name"].GetString()]
        self.name = parameters["name"].GetString()

        sensitivity_variable_names = parameters["sensitivity_variable_names"].GetStringArray()

        allowed_sensitivity_variables = [
            "SHAPE_SENSITIVITY",
            "DENSITY_SENSITIVITY",
            "THICKNESS_SENSITIVITY",
            "CROSS_AREA_SENSITIVITY"
        ]

        self.sensitivity_variables = []
        for sensitivity_variable_name in sensitivity_variable_names:
            if sensitivity_variable_name not in allowed_sensitivity_variables:
                raise RuntimeError(f"{sensitivity_variable_name} used in {self.name} is not supported. Followings are supported \"sensitivity_variable_names\":\n\t" + "\n\t".join(allowed_sensitivity_variables))

            self.sensitivity_variables.append(Kratos.KratosGlobals.GetVariable(sensitivity_variable_name))

    def Check(self, _: dict):
        data_communicator = self.model_part.GetCommunicator().GetDataCommunicator()
        if not KratosOA.OptimizationVariableUtils.IsVariableExistsInAllContainerProperties(self.model_part.Elements, Kratos.DENSITY, data_communicator):
            raise RuntimeError(f"Some elements' properties in {self.model_part.FullName()} does not have DENSITY variable.")

        if KratosOA.OptimizationVariableUtils.IsVariableExistsInAtLeastOneContainerProperties(self.model_part.Elements, Kratos.THICKNESS, data_communicator) and \
           KratosOA.OptimizationVariableUtils.IsVariableExistsInAtLeastOneContainerProperties(self.model_part.Elements, Kratos.CROSS_AREA, data_communicator):
           raise RuntimeError(f"{self.model_part.FullName()} has elements consisting THICKNESS and CROSS_AREA. Please break down this response to SumResponseFunction where each sub response function only has elements with either THICKNESS or CROSS_AREA.")

        if not KratosOA.OptimizationVariableUtils.AreAllEntitiesOfSameGeometryType(self.model_part.Elements, data_communicator):
            raise RuntimeError(f"{self.model_part.FullName()} has elements with different geometry types. Please break down this response to SumResponseFunction where each sub response function only has elements with one geometry type.")

    def Initialize(self, _: dict):
        pass

    def InitializeIteration(self, _: dict):
        pass

    def FinalizeIteration(self, _: dict):
        pass

    def Finalize(self, _: dict):
        pass

    def CalculateValue(self, _: dict):
        return KratosOA.MassResponseUtilities.CalculateMass(self.model_part)

    def CalculateSensitivity(self, _: dict):
        response_sensitivities = []
        for sensitivity_variable in self.sensitivity_variables:
            if sensitivity_variable == Kratos.SHAPE_SENSITIVITY:
                KratosOA.MassResponseUtilities.CalculateMassShapeSensitivity(self.model_part, sensitivity_variable)
                sensitivity_container_enum = ContainerEnum.NODES
            elif sensitivity_variable == KratosOA.DENSITY_SENSITIVITY:
                KratosOA.MassResponseUtilities.CalculateMassDensitySensitivity(self.model_part, sensitivity_variable)
                sensitivity_container_enum = ContainerEnum.ELEMENT_PROPERTIES
            elif sensitivity_variable == KratosOA.THICKNESS_SENSITIVITY:
                KratosOA.MassResponseUtilities.CalculateMassThicknessSensitivity(self.model_part, sensitivity_variable)
                sensitivity_container_enum = ContainerEnum.ELEMENT_PROPERTIES
            elif sensitivity_variable == KratosOA.CROSS_AREA_SENSITIVITY:
                KratosOA.MassResponseUtilities.CalculateMassCrossAreaSensitivity(self.model_part, sensitivity_variable)
                sensitivity_container_enum = ContainerEnum.ELEMENT_PROPERTIES
            else:
                RuntimeError(f"Unsupported sensitivity {sensitivity_variable.Name} requested in response with name {self.name}.")

            response_sensitivities.append(ResponseSentivity(sensitivity_variable, sensitivity_container_enum, self.model_part))

        return response_sensitivities

    def GetResponseFunctionName(self):
        return self.name

