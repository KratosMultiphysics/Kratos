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
            "name"                 : "",
            "echo_level"           : 0,
            "module"               : "KratosMultiphysics.OptimizationApplication.controls",
            "type"                 : "PleaseProvideClassName",
            "sensitivity_modifiers": [],
            "update_modifiers"     : [],
            "settings"             : {}
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

        modifier_defaults = Kratos.Parameters("""{
            "module"  : "KratosMultiphysics.OptimizationApplication.modifiers",
            "type"    : "PleaseProvideClassName",
            "settings": {}
        }""")

        self.__list_of_sensitivity_modifiers = []
        for modifier_settings in self.parameters["sensitivity_modifiers"]:
            modifier_settings.ValidateAndAssignDefaults(modifier_defaults)
            self.__list_of_sensitivity_modifiers.append(RetrieveObject(self.model, modifier_settings, optimization_info, Modifier))

        self.__list_of_update_modifiers = []
        for modifier_settings in self.parameters["update_modifiers"]:
            modifier_settings.ValidateAndAssignDefaults(modifier_defaults)
            self.__list_of_update_modifiers.append(RetrieveObject(self.model, modifier_settings, optimization_info, Modifier))

    def Initialize(self):
        self.control.Initialize()
        for modifier in self.__list_of_sensitivity_modifiers:
            modifier.Initialize()

        for modifier in self.__list_of_update_modifiers:
            modifier.Initialize()

    def InitializeSolutionStep(self):
        self.control.InitializeSolutionStep()

        for modifier in self.__list_of_sensitivity_modifiers:
            modifier.InitializeSolutionStep()

        for modifier in self.__list_of_update_modifiers:
            modifier.InitializeSolutionStep()

    def FinalizeSolutionStep(self):
        self.control.FinalizeSolutionStep()

        for modifier in self.__list_of_sensitivity_modifiers:
            modifier.FinalizeSolutionStep()

        for modifier in self.__list_of_update_modifiers:
            modifier.FinalizeSolutionStep()

    def Finalize(self):
        self.control.Finalize()

        for modifier in self.__list_of_sensitivity_modifiers:
            modifier.Finalize()

        for modifier in self.__list_of_update_modifiers:
            modifier.Finalize()

    def GetName(self):
        return self.name

    def GetControl(self) -> Control:
        return self.control

    def ModifySensitivities(self, sensitivities: Kratos.Vector) -> Kratos.Vector:
        resultant_sensitivities = Kratos.Vector(sensitivities)
        for modifier in self.__list_of_sensitivity_modifiers:
            resultant_sensitivities = modifier.ModifySensitivities(resultant_sensitivities, self.control.GetModelPart(), self.control.GetContainerType())

        return resultant_sensitivities

    def ModifyControlUpdates(self, controls: Kratos.Vector) -> Kratos.Vector:
        resultant_controls = Kratos.Vector(controls)
        for modifier in self.__list_of_update_modifiers:
            resultant_controls = modifier.ModifySensitivities(resultant_controls, self.control.GetModelPart(), self.control.GetContainerType())

        return resultant_controls

