import KratosMultiphysics as Kratos
from KratosMultiphysics.OptimizationApplication.optimization_routine import OptimizationRoutine
from KratosMultiphysics.OptimizationApplication.optimization_info import OptimizationInfo
from KratosMultiphysics.OptimizationApplication.utilities.helper_utils import Factory
from KratosMultiphysics.OptimizationApplication.utilities.helper_utils import ContainerVariableDataHolder
from KratosMultiphysics.OptimizationApplication.controls.control import Control
from KratosMultiphysics.OptimizationApplication.modifiers.modifier import Modifier

class ControlWrapper(OptimizationRoutine):
    def __init__(self, model: Kratos.Model, parameters: Kratos.Parameters, optimization_info: OptimizationInfo):
        self.model = model
        self.parameters = parameters
        self.optimization_info = optimization_info
        self.__control_updates = {}

        default_parameters = Kratos.Parameters("""{
            "name"          : "",
            "echo_level"    : 0,
            "module"        : "KratosMultiphysics.OptimizationApplication.controls",
            "type"          : "PleaseProvideClassName",
            "modifiers_list": [],
            "settings"      : {}
        }""")

        self.parameters.ValidateAndAssignDefaults(default_parameters)

        self.name = self.parameters["name"].GetString()
        if self.name == "":
            raise RuntimeError(f"Controls should be given a non-empty name. Followings are the corresponding control settings:\n{self.parameters}")

        # create the response function
        self.control: Control = Factory(self.parameters["module"].GetString(), self.parameters["type"].GetString(), self.model, self.parameters["settings"], optimization_info, Control)

        self.__list_of_modifiers = []
        for modifier_settings in self.parameters["modifiers_list"]:
            self.__list_of_modifiers.append(Factory(self.model, modifier_settings, optimization_info, Modifier))

    def Initialize(self):
        self.control.Initialize()
        for modifier in self.__list_of_modifiers:
            modifier.Initialize()

    def InitializeSolutionStep(self):
        if self.optimization_info["step"] > 1:
            control_update: ContainerVariableDataHolder
            for control_update in self.__control_updates.values():
                current_control_values = self.control.GetCurrentControls(control_update.GetModelPart())
                new_control_values = control_update.GetValues() + current_control_values.GetValues()
                modified_control_values = self.ModifyControl(ContainerVariableDataHolder(new_control_values, current_control_values.GetVariable(), current_control_values.GetModelPart(), current_control_values.GetContainerType()))
                self.control.UpdateControls(modified_control_values)

        self.__control_updates = {}

        self.control.InitializeSolutionStep()

        for modifier in self.__list_of_modifiers:
            modifier.InitializeSolutionStep()

    def FinalizeSolutionStep(self):
        self.control.FinalizeSolutionStep()

        for modifier in self.__list_of_modifiers:
            modifier.FinalizeSolutionStep()

    def Finalize(self):
        self.control.Finalize()

        for modifier in self.__list_of_modifiers:
            modifier.Finalize()

    def GetName(self):
        return self.name

    def GetControl(self) -> Control:
        return self.control

    def ModifyControl(self, control_values: ContainerVariableDataHolder):
        resultant_control_values = ContainerVariableDataHolder(Kratos.Vector(control_values.GetValues()), control_values.GetVariable(), control_values.GetModelPart(), control_values.GetContainerType())
        for modifier in self.__list_of_modifiers:
            resultant_control_values = modifier.Control(resultant_control_values, self.control.GetModelPart(), self.control.GetContainerType())

        return resultant_control_values

    def ModifySensitivities(self, sensitivities: ContainerVariableDataHolder) -> ContainerVariableDataHolder:
        resultant_sensitivities = ContainerVariableDataHolder(Kratos.Vector(sensitivities.GetValues()), sensitivities.GetVariable(), sensitivities.GetModelPart(), sensitivities.GetContainerType())
        for modifier in self.__list_of_modifiers:
            resultant_sensitivities = modifier.ModifySensitivities(resultant_sensitivities, self.control.GetModelPart(), self.control.GetContainerType())

        return resultant_sensitivities

    def ModifyControlUpdates(self, controls: ContainerVariableDataHolder) -> ContainerVariableDataHolder:
        resultant_controls = ContainerVariableDataHolder(Kratos.Vector(controls.GetValues()), controls.GetVariable(), controls.GetModelPart(), controls.GetContainerType())
        for modifier in self.__list_of_modifiers:
            resultant_controls = modifier.ModifyControlUpdates(resultant_controls, self.control.GetModelPart(), self.control.GetContainerType())

        return resultant_controls

    def SetControlUpdate(self, control_update: ContainerVariableDataHolder):
        # check whether it is a valid control update
        control_model_part: Kratos.ModelPart
        for control_model_part in self.control.GetModelParts():
            is_valid_update = True

            is_valid_update = is_valid_update and (control_model_part.FullName() in [model_part.FullName() for model_part in self.control.GetModelParts()])
            is_valid_update = is_valid_update and (control_update.GetContainerType() == self.control.GetContainerType())
            is_valid_update = is_valid_update and (control_update.GetVariable() == self.control.GetControlUpdateVariable())

            if is_valid_update:
                break

        if is_valid_update:
            control_update_key = (control_model_part.FullName(), control_update.GetContainerType(), control_update.GetVariable())
            if control_update_key in self.__control_updates.keys():
                raise RuntimeError(f"Trying to overwrite control update for {self.GetName()} with {str(control_update)}")

            # set the control update
            self.__control_updates[control_update_key] = control_update
        else:
            raise RuntimeError(f"Trying to set an invalid control update for {self.GetName} with {str(control_update)}")

