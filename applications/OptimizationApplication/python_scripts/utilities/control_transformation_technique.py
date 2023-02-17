import KratosMultiphysics as Kratos
import KratosMultiphysics.OptimizationApplication as KratosOA
from KratosMultiphysics.OptimizationApplication.optimization_info import OptimizationInfo
from KratosMultiphysics.OptimizationApplication.controls.control import Control
from KratosMultiphysics.OptimizationApplication.transformation_techniques.transformation_technique import TransformationTechnique
from KratosMultiphysics.OptimizationApplication.utilities.helper_utilities import OptimizationProcessFactory
from KratosMultiphysics.OptimizationApplication.utilities.helper_utilities import CallOnAll
from KratosMultiphysics.OptimizationApplication.utilities.helper_utilities import WriteCollectiveVariableDataHolderToOptmizationInfo

class ControlTransformationTechnique(Kratos.Process):
    def __init__(self, model: Kratos.Model, parameters: Kratos.Parameters, optimization_info: OptimizationInfo):
        super().__init__()

        default_parameters = Kratos.Parameters("""{
            "name"                     : "",
            "module"                   : "KratosMultiphysics.OptimizationApplication.controls",
            "type"                     : "PleaseProvideClassName",
            "transformation_techniques": [],
            "settings"                 : {}
        }""")

        parameters.ValidateAndAssignDefaults(default_parameters)

        self.__name = parameters["name"].GetString()
        self.__control: Control = OptimizationProcessFactory(parameters["module"].GetString(), parameters["type"].GetString(), model, parameters["settings"], optimization_info, Control)
        self.__control.SetName(self.GetName())
        self.__transformation_techniques: 'list[TransformationTechnique]' = []
        self.__optimization_info = optimization_info

        default_transformation_settings = Kratos.Parameters("""{
            "module"  : "KratosMultiphysics.OptimizationApplication.transformation_techniques",
            "type"    : "PleaseProvideClassName",
            "settings": {}
        }""")
        for transformation_technique_settings in parameters["transformation_techniques"]:
            transformation_technique_settings.ValidateAndAssignDefaults(default_transformation_settings)
            self.__transformation_techniques.append(OptimizationProcessFactory(transformation_technique_settings["module"].GetString(), transformation_technique_settings["type"].GetString(), model, transformation_technique_settings["settings"], optimization_info, TransformationTechnique))

    def GetName(self) -> str:
        return self.__name

    def GetControl(self):
        return self.__control

    def GetTransformationTechniques(self) -> 'list[TransformationTechnique]':
        return self.__transformation_techniques

    def TransformSensitivity(self, collective_container_variable_data_holder: KratosOA.CollectiveVariableDataHolder):
        for container in collective_container_variable_data_holder.GetVariableDataHolders():
            CallOnAll(self.__transformation_techniques, TransformationTechnique.TransformSensitivity, container)

    def TransformUpdate(self, collective_container_variable_data_holder: KratosOA.CollectiveVariableDataHolder):
        for container in collective_container_variable_data_holder.GetVariableDataHolders():
            CallOnAll(self.__transformation_techniques, TransformationTechnique.TransformUpdate, container)

    def SetControlUpdate(self, collective_container_variable_data_holder: KratosOA.CollectiveVariableDataHolder):
        WriteCollectiveVariableDataHolderToOptmizationInfo(
            self.__optimization_info,
            collective_container_variable_data_holder.Clone(),
            f"problem_data/control_data/<model_part_name>/{self.GetName()}/update")

        self.__optimization_info.SetValue(f"problem_data/control_data/{self.GetName()}/collective_update", collective_container_variable_data_holder)

    def ApplyControlUpdate(self):
        key = f"problem_data/control_data/{self.GetName()}/collective_update"
        if self.__optimization_info.HasValue(key):
            self.GetControl().UpdateControl(self.__optimization_info.GetValue(key))
        else:
            raise RuntimeError(f"Control update not set in {self.GetName()} controller technique.")