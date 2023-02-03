from typing import Union

import KratosMultiphysics as Kratos
import KratosMultiphysics.OptimizationApplication as KratosOA

from KratosMultiphysics.OptimizationApplication.optimization_info import OptimizationInfo
from KratosMultiphysics.OptimizationApplication.controls.control import Control
from KratosMultiphysics.OptimizationApplication.transformation_techniques.transformation_technique import TransformationTechnique
from KratosMultiphysics.OptimizationApplication.utilities.helper_utils import Factory
from KratosMultiphysics.OptimizationApplication.utilities.helper_utils import CallOnAll
from KratosMultiphysics.OptimizationApplication.utilities.helper_utils import ContainerVariableDataHolderUnion

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

    def TransformSensitivity(self, container_variable_data_holder: ContainerVariableDataHolderUnion):
        CallOnAll(self.__transformation_techniques, TransformationTechnique.TransformSensitivity, container_variable_data_holder)

    def TransformUpdate(self, container_variable_data_holder: ContainerVariableDataHolderUnion):
        CallOnAll(self.__transformation_techniques, TransformationTechnique.TransformUpdate, container_variable_data_holder)

    def SetControlUpdate(self, container_variable_data_holder: ContainerVariableDataHolderUnion):
        # check whether container data is valid for the updates
        if container_variable_data_holder.GetModelPart() not in self.GetControl().GetModelParts():
            raise RuntimeError(f"The control update for {container_variable_data_holder.GetModelPart().FullName()} is not one of the controlled model parts in control with name \"{self.GetName()}\". Followings are allowd model parts: \n\t" + "\n\t".join([v.FullName() for v in self.GetControl().GetModelParts()]))

        if not isinstance(container_variable_data_holder, self.GetControl().CreateContainerVariableDataHolder(container_variable_data_holder.GetModelPart()).__class__):
            raise RuntimeError(f"The control type {container_variable_data_holder.__class__.__name__} mismatch with the required control type {self.GetControl().CreateContainerVariableDataHolder(container_variable_data_holder.GetModelPart()).__class__.__name__} in \"{self.GetName()}\" control.")

        self.__control_updates[container_variable_data_holder.GetModelPart()] = container_variable_data_holder.Clone()

    def ApplyControlUpdate(self):
        for model_part in self.GetControl().GetModelParts():
            if model_part not in self.__control_updates.keys():
                raise RuntimeError(f"Control update for {model_part.FullName()} not set in {self.GetName()} controller technique.")

            self.GetControl().UpdateControl(self.__control_updates[model_part])