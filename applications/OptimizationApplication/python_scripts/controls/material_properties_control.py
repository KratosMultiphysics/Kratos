import KratosMultiphysics as Kratos
import KratosMultiphysics.OptimizationApplication as KratosOA
from KratosMultiphysics.OptimizationApplication.controls.control import Control
from KratosMultiphysics.OptimizationApplication.optimization_info import OptimizationInfo
from KratosMultiphysics.OptimizationApplication.utilities.container_data import ContainerData

class MaterialPropertiesControl(Control):
    def __init__(self, model: Kratos.Model, parameters: Kratos.Parameters, optimization_info: OptimizationInfo):
        default_settings = Kratos.Parameters("""{
            "model_part_names"            : [""],
            "control_variable_name"       : "",
            "control_update_variable_name": "SCALAR_CONTROL_UPDATE"
        }""")
        parameters.ValidateAndAssignDefaults(default_settings)

        self.optimization_info = optimization_info
        self.model_parts = [model[model_part_name] for model_part_name in parameters["model_part_names"].GetStringArray()]

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
        # creates element specific properties
        if not self.optimization_info.HasSolutionStepDataKey("model_parts_with_element_specific_properties"):
            self.optimization_info["model_parts_with_element_specific_properties"] = []

        for model_part in self.model_parts:
            if not f"{model_part.FullName()}.Elements" in self.optimization_info["model_parts_with_element_specific_properties"]:
                KratosOA.OptimizationUtils.CreateEntitySpecificPropertiesForContainer(model_part, model_part.Elements)
                self.optimization_info["model_parts_with_element_specific_properties"].append(f"{model_part.FullName()}.Elements")

    def UpdateControl(self, control_values: ContainerData):
        current_values_container = ContainerData(control_values.GetModelPart(), control_values.GetContainerTpe())
        current_values_container.ReadDataFromContainerVariable(self.control_variable)
        new_values_container = current_values_container + control_values
        new_values_container.AssignDataToContainer(self.control_variable)

    def GetModelParts(self) -> 'list[Kratos.ModelPart]':
        return self.model_parts

    def GetContainerType(self) -> ContainerData.ContainerEnum:
        return ContainerData.ContainerEnum.ELEMENT_PROPERTIES

    def GetControlSensitivityVariable(self) -> any:
        return self.sensitivity_variable

    def GetControlUpdateVariable(self) -> any:
        return self.control_update_variable

