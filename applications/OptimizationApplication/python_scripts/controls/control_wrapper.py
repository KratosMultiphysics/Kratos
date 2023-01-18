import KratosMultiphysics as Kratos
from KratosMultiphysics.OptimizationApplication.optimization_routine import OptimizationRoutine
from KratosMultiphysics.OptimizationApplication.optimization_info import OptimizationInfo
from KratosMultiphysics.OptimizationApplication.utilities.helper_utils import Factory
from KratosMultiphysics.OptimizationApplication.utilities.container_data import ContainerData
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
        # update the controls is the step is > 1
        if self.optimization_info["step"] > 1:
            control_update: ContainerData
            for control_update in self.__control_updates.values():
                current_control_values = self.control.GetCurrentControlContainerData(control_update.GetModelPart())
                new_control_values = current_control_values + control_update
                self.ModifyControl(new_control_values)
                self.control.UpdateControls(new_control_values)

        # reset control updates
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

    def ModifyControl(self, controls: ContainerData):
        for modifier in self.__list_of_modifiers:
             modifier.Control(controls)

    def ModifySensitivities(self, sensitivities: ContainerData):
        for modifier in self.__list_of_modifiers:
            modifier.ModifySensitivities(sensitivities)

    def ModifyControlUpdates(self, control_update: ContainerData):
        for modifier in self.__list_of_modifiers:
            modifier.ModifyControlUpdates(control_update)

    def SetControlUpdate(self, control_update: ContainerData):
        control_model_part: Kratos.ModelPart
        for control_model_part in self.control.GetModelParts():
            if control_update.IsSameContainer(ContainerData(control_model_part, self.GetControl().GetContainerType())):
                self.__control_updates[control_model_part] = control_update
                return

        raise RuntimeError(f"Trying to set an invalid control update for \"{self.GetName()}\" with {str(control_update)}")

