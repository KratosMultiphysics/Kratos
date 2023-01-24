import KratosMultiphysics as Kratos

from KratosMultiphysics.OptimizationApplication.optimization_info import OptimizationInfo
from KratosMultiphysics.OptimizationApplication.controls.control import Control
from KratosMultiphysics.OptimizationApplication.transformation_techniques.transformation_technique import TransformationTechnique
from KratosMultiphysics.OptimizationApplication.utilities.container_data import ContainerData
from KratosMultiphysics.OptimizationApplication.utilities.helper_utils import Factory
from KratosMultiphysics.OptimizationApplication.utilities.helper_utils import CallOnAll

class ControlTransformationTechnique:
    def __init__(self, model: Kratos.Model, parameters: Kratos.Parameters, optimization_info: OptimizationInfo):
        default_parameters = Kratos.Parameters("""{
            "name"                     : "",
            "module"                   : "KratosMultiphysics.OptimizationApplication.controls",
            "type"                     : "PleaseProvideClassName",
            "transformation_techniques": [],
            "settings"                 : {}
        }""")

        parameters.ValidateAndAssignDefaults(default_parameters)

        self.__name = parameters["name"].GetString()
        self.__control: Control = Factory(parameters["module"].GetString(), parameters["type"].GetString(), model, parameters["settings"], optimization_info, Control)
        self.__transformation_techniques: 'list[TransformationTechnique]' = []
        self.__control_updates = {}

        default_transformation_settings = Kratos.Parameters("""{
            "module"  : "KratosMultiphysics.OptimizationApplication.transformation_techniques",
            "type"    : "PleaseProvideClassName",
            "settings": {}
        }""")
        for transformation_technique_settings in parameters["transformation_techniques"]:
            transformation_technique_settings.ValidateAndAssignDefaults(default_transformation_settings)
            self.__transformation_techniques.append(Factory(transformation_technique_settings["module"].GetString(), transformation_technique_settings["type"].GetString(), model, transformation_technique_settings["settings"], optimization_info, TransformationTechnique))

    def GetName(self) -> str:
        return self.__name

    def GetControl(self):
        return self.__control

    def Initialize(self):
        self.__control.Initialize()
        CallOnAll(self.__transformation_techniques, TransformationTechnique.Initialize)

    def InitializeSolutionStep(self):
        self.__control_updates = {}
        self.__control.InitializeSolutionStep()
        CallOnAll(self.__transformation_techniques, TransformationTechnique.InitializeSolutionStep)

    def FinalizeSolutionStep(self):
        self.__control.FinalizeSolutionStep()
        CallOnAll(self.__transformation_techniques, TransformationTechnique.FinalizeSolutionStep)

    def Finalize(self):
        self.__control.Finalize()
        CallOnAll(self.__transformation_techniques, TransformationTechnique.Finalize)

    def TransformSensitivity(self, container_data: ContainerData):
        CallOnAll(self.__transformation_techniques, TransformationTechnique.Finalize, container_data)

    def TransformUpdate(self, container_data: ContainerData):
        CallOnAll(self.__transformation_techniques, TransformationTechnique.Finalize, container_data)

    def SetControlUpdate(self, container_data: ContainerData):
        # check whether container data is valid for the updates
        if container_data.GetModelPart() not in self.GetControl().GetModelParts():
            raise RuntimeError(f"The control update for {container_data.GetModelPart().FullName()} is not one of the controlled model parts in control with name \"{self.GetName()}\". Followings are allowd model parts: \n\t" + "\n\t".join([v.FullName() for v in self.GetControl().GetModelParts()]))

        if container_data.GetContainerTpe() != self.GetControl().GetContainerType():
            raise RuntimeError(f"The control type {container_data.GetContainerTpe().name} mismatch with the required control type {self.GetControl().GetContainerType().name} in \"{self.GetName()}\" control.")

        self.__control_updates[container_data.GetModelPart()] = container_data

    def ApplyControlUpdate(self):
        for model_part in self.GetControl().GetModelParts():
            if model_part not in self.__control_updates.keys():
                raise RuntimeError(f"Control update for {model_part.FullName()} not set in {self.GetName()} controller technique.")

            self.GetControl().UpdateControl(self.__control_updates[model_part])