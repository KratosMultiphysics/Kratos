import KratosMultiphysics as Kratos
import KratosMultiphysics.OptimizationApplication as KratosOA
from KratosMultiphysics.OptimizationApplication.controls.control import Control
from KratosMultiphysics.OptimizationApplication.optimization_info import OptimizationInfo
from KratosMultiphysics.OptimizationApplication.utilities.helper_utils import ContainerEnum

class PropertiesControl(Control):
    def __init__(self, model: Kratos.Model, parameters: Kratos.Parameters, optimization_info: OptimizationInfo):
        super().__init__(model, parameters, optimization_info)

        default_settings = Kratos.Parameters("""{
            "model_part_name"             : "",
            "control_variable_name"       : "",
            "control_update_variable_name": "SCALAR_CONTROL_UPDATE"
        }""")
        parameters.ValidateAndAssignDefaults(default_settings)

        self.model_part = self.model[parameters["model_part_name"].GetString()]

        control_variable_name = parameters["control_variable_name"].GetString()
        control_variable_type = Kratos.KratosGlobals.GetVariableType(control_variable_name)
        if control_variable_type != "Double":
            raise RuntimeError(f"{control_variable_name} with {control_variable_type} type is not supported. Only supports double variables")

        sensitivity_variable_name = control_variable_name + "_SENSITIVITY"
        if not Kratos.KratosGlobals.HasVariable(sensitivity_variable_name):
            raise RuntimeError(f"Sensitivity variable {sensitivity_variable_name} for control variable {control_variable_name} not found.")
        elif Kratos.KratosGlobals.GetVariableType(sensitivity_variable_name) != "Double":
            raise RuntimeError(f"Sensitivity variable {sensitivity_variable_name} type is not double.")

        control_update_variable_name = parameters["control_update_variable_name"].GetString()
        control_update_variable_type = Kratos.KratosGlobals.GetVariableType(control_update_variable_name)
        if control_update_variable_type != "Double":
            raise RuntimeError(f"{control_update_variable_name} with {control_update_variable_type} type is not supported. Only supports double variables")

        self.control_variable = Kratos.KratosGlobals.GetVariable(control_variable_name)
        self.sensitivity_variable = Kratos.KratosGlobals.GetVariable(sensitivity_variable_name)
        self.control_update_variable = Kratos.KratosGlobals.GetVariable(control_update_variable_name)

    def Initialize(self):
        super().Initialize()

        # creates element specific properties
        if not self.optimization_info.HasSolutionStepDataKey("model_parts_with_element_specific_properties"):
            self.optimization_info["model_parts_with_element_specific_properties"] = []

        if not self.model_part.FullName() in self.optimization_info["model_parts_with_element_specific_properties"]:
            KratosOA.OptimizationUtils.CreateEntitySpecificPropertiesForContainer(self.model_part, self.model_part.Elements)
            self.optimization_info["model_parts_with_element_specific_properties"].append(self.model_part.FullName())

    def UpdateControls(self, control_values: Kratos.Vector):
        KratosOA.OptimizationUtils.AssignVectorToContainerProperties(self.model_part.Elements, self.control_variable, control_values)
        KratosOA.OptimizationUtils.AssignVectorToContainer(self.model_part.Elements, self.model_part.ProcessInfo[Kratos.DOMAIN_SIZE], self.control_variable, control_values)

    def GetNewControlValuesVector(self) -> Kratos.Vector:
        current_property_values = Kratos.Vector()
        KratosOA.OptimizationUtils.GetContainerPropertiesVariableToVector(self.model_part.Elements, self.control_variable, current_property_values)
        return current_property_values + self.GetControlUpdatesVector()

    def GetModelPart(self) -> Kratos.ModelPart:
        return self.model_part

    def GetContainerType(self) -> ContainerEnum:
        return ContainerEnum.ELEMENT_PROPERTIES

    def GetControlSensitivityVariable(self) -> any:
        return self.sensitivity_variable

    def GetControlUpdateVariable(self) -> any:
        return self.control_update_variable

