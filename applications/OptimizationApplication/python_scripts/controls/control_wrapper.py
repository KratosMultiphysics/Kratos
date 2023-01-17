import KratosMultiphysics as Kratos
from KratosMultiphysics.OptimizationApplication.optimization_routine import OptimizationRoutine
from KratosMultiphysics.OptimizationApplication.optimization_info import OptimizationInfo
from KratosMultiphysics.OptimizationApplication.utilities.helper_utils import RetrieveObject
from KratosMultiphysics.OptimizationApplication.controls.control import Control
from KratosMultiphysics.OptimizationApplication.modifiers.modifier import Modifier

class ControlWrapper(OptimizationRoutine):
    def __init__(self, model: Kratos.Model, parameters: Kratos.Parameters, optimization_info: OptimizationInfo):
        self.model = model
        self.parameters = parameters
        self.optimization_info = optimization_info

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
        control_settings = Kratos.Parameters("""{}""")
        control_settings.AddString("module", self.parameters["module"].GetString())
        control_settings.AddString("type", self.parameters["type"].GetString())
        control_settings.AddValue("settings", self.parameters["settings"])
        self.control: Control = RetrieveObject(self.model, control_settings, optimization_info, Control)

        self.__list_of_modifiers = []
        for modifier_settings in self.parameters["modifiers_list"]:
            self.__list_of_modifiers.append(RetrieveObject(self.model, modifier_settings, optimization_info, Modifier))

    def Initialize(self):
        self.control.Initialize()
        for modifier in self.__list_of_modifiers:
            modifier.Initialize()

    def InitializeSolutionStep(self):
        if self.optimization_info["step"] > 1:
            control_values = self.control.GetNewControlValuesVector()
            control_values = self.ModifyControl(control_values)
            self.control.UpdateControls(control_values)

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

    def ModifyControl(self, control_values: Kratos.Vector):
        resultant_control_values = Kratos.Vector(control_values)
        for modifier in self.__list_of_modifiers:
            resultant_control_values = modifier.Control(resultant_control_values, self.control.GetModelPart(), self.control.GetContainerType())

        return resultant_control_values

    def ModifySensitivities(self, sensitivities: Kratos.Vector) -> Kratos.Vector:
        resultant_sensitivities = Kratos.Vector(sensitivities)
        for modifier in self.__list_of_modifiers:
            resultant_sensitivities = modifier.ModifySensitivities(resultant_sensitivities, self.control.GetModelPart(), self.control.GetContainerType())

        return resultant_sensitivities

    def ModifyControlUpdates(self, controls: Kratos.Vector) -> Kratos.Vector:
        resultant_controls = Kratos.Vector(controls)
        for modifier in self.__list_of_modifiers:
            resultant_controls = modifier.ModifyControlUpdates(resultant_controls, self.control.GetModelPart(), self.control.GetContainerType())

        return resultant_controls

