import KratosMultiphysics as Kratos
import KratosMultiphysics.OptimizationApplication as KratosOA
from KratosMultiphysics.OptimizationApplication.controls.control import Control
from KratosMultiphysics.OptimizationApplication.optimization_info import OptimizationInfo
from KratosMultiphysics.OptimizationApplication.utilities.response_function_implementor import ResponseFunctionImplementor

class MaterialPropertiesControl(Control):
    def __init__(self, model: Kratos.Model, parameters: Kratos.Parameters, optimization_info: OptimizationInfo):
        super().__init__()

        default_settings = Kratos.Parameters("""{
            "model_part_names"            : [""],
            "control_variable_name"       : ""
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

        self.control_variable = Kratos.KratosGlobals.GetVariable(control_variable_name)
        self.sensitivity_variable = Kratos.KratosGlobals.GetVariable(sensitivity_variable_name)

    def ExecuteInitialize(self):
        # creates element specific properties
        if "model_parts_with_element_specific_properties" not in self.optimization_info.GetSolutionStepData(0).keys():
            self.optimization_info["model_parts_with_element_specific_properties"] = []

        for model_part in self.model_parts:
            if not f"{model_part.FullName()}.Elements" in self.optimization_info["model_parts_with_element_specific_properties"]:
                KratosOA.OptimizationUtils.CreateEntitySpecificPropertiesForContainer(model_part, model_part.Elements)
                self.optimization_info["model_parts_with_element_specific_properties"].append(f"{model_part.FullName()}.Elements")

    def CalculateSensitivity(self, response_function: ResponseFunctionImplementor, output_sensitivities: KratosOA.CollectiveVariableDataHolder):
        # clear the container
        output_sensitivities.ClearVariableDataHolders()

        for model_part in self.model_parts:
            d_j_d_phi = KratosOA.ElementPropertiesContainerVariableDataHolder(model_part)
            response_function.GetStandardizedSensitivity(self.sensitivity_variable, d_j_d_phi)

            output_sensitivities.AddVariableDataHolder(d_j_d_phi)

    def UpdateControl(self, update: KratosOA.CollectiveVariableDataHolder):
        for current_update in update.GetVariableDataHolders():
            model_part: Kratos.ModelPart = current_update.GetModelPart()

            if model_part not in self.model_parts:
                raise RuntimeError(f"Trying to update {model_part.FullName()} which is not controlled by {self.GetName()}.")

            current_values_container = KratosOA.ElementPropertiesContainerVariableDataHolder(model_part)
            current_values_container.ReadDataFromContainerVariable(self.control_variable)
            new_values_container: KratosOA.ElementPropertiesContainerVariableDataHolder = current_values_container + current_update
            new_values_container.AssignDataToContainerVariable(self.control_variable)

