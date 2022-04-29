import KratosMultiphysics as KM
from KratosMultiphysics import Parameters, Logger
import KratosMultiphysics.OptimizationApplication as KOA
from KratosMultiphysics.OptimizationApplication.responses.base_response import BaseResponseFunction
import KratosMultiphysics.StructuralMechanicsApplication as KSM

import time as timer
import numpy as np

class LinearResponseFunction(BaseResponseFunction):

    def __init__(self,response_name, response_settings,model):

        self.type = "linear"
        self.variable = "LINEAR"
        super().__init__(response_name, response_settings, model)

        if not self.response_settings.Has("gradient_settings"):
            self.gradient_settings = KM.Parameters()
            self.gradient_settings.AddString("gradient_mode","semi_analytic")
            self.gradient_settings.AddDouble("step_size",1e-6)
        else:
            self.gradient_settings = self.response_settings["gradient_settings"]     

        self.supported_control_types = ["shape"]
        self.gradients_variables = {"shape":"D_LINEAR_D_X"}

        if len(self.evaluated_model_parts) != 1:
            raise RuntimeError("LinearResponseFunction: 'evaluated_objects' of response '{}' must have only one entry !".format(self.name)) 

        for control_type in self.control_types:
            if not control_type in self.supported_control_types:
                raise RuntimeError("LinearResponseFunction: type {} in 'control_types' of response '{}' is not supported, supported types are {}  !".format(control_type,self.name,self.supported_control_types)) 

        
        root_model_part_name = self.evaluated_model_parts[0].split(".")[0]
        for evaluated_model_part in self.evaluated_model_parts:
            if evaluated_model_part.split(".")[0] != root_model_part_name:
                raise RuntimeError("LinearResponseFunction: evaluated_model_parts of mass response must have the same root model part !")

        self.root_model_part = self.model.GetModelPart(root_model_part_name)
        # add vars
        for control_type in self.control_types:
            if control_type == "shape":
                self.root_model_part.AddNodalSolutionStepVariable(KM.KratosGlobals.GetVariable(self.gradients_variables[control_type]))

    def GetVariableName(self):
        return  self.variable

    def GetGradientsVariablesName(self):
        return self.gradients_variables

    def GetGradientVariableNameForType(self,control_type, raise_error=True):
        if raise_error:
            if not control_type in self.supported_control_types:
                raise RuntimeError("LinearResponseFunction: type {} in 'control_types' of response '{}' is not supported, supported types are {}  !".format(control_type,self.name,self.supported_control_types)) 

        return self.gradients_variables[control_type]

    def Initialize(self):
        super().Initialize()

    def CalculateValue(self):
        Logger.PrintInfo("LinearResponseFunction:CalculateValue: Starting value calculation for response ", self.name)
        startTime = timer.time()
        self.value = 1
        Logger.PrintInfo("LinearResponseFunction:CalculateValue: Time needed for calculating value ",round(timer.time() - startTime,2),"s")        
        return self.value

    def CalculateGradientsForTypesAndObjects(self,control_types,controlled_objects,raise_error=True):

        if raise_error:
            for itr in range(len(controlled_objects)):
                controlled_object = controlled_objects[itr]
                control_type = control_types[itr]
                found = False
                for itr_2 in range(len(self.controlled_model_parts)):
                    controlled_model_part = self.controlled_model_parts[itr_2]
                    controlled_type = self.control_types[itr_2]
                    if controlled_type==control_type and controlled_model_part==controlled_object:
                        found = True
                        break
                if not found:
                    raise RuntimeError("LinearResponseFunction:CalculateGradientsForTypesAndObjects: control type {} of control object {} is not in the control_types of response {}".format(control_types[itr],controlled_object,self.name))

        Logger.PrintInfo("LinearResponseFunction", "Starting ", control_types," gradients calculation of response ", self.name," for ",controlled_objects)
        startTime = timer.time()

        controlle_object_index = 0
        for controlle_object in controlled_objects:
            model_part = self.model.GetModelPart(controlle_object)
            control_type = control_types[controlle_object_index]           
            for node in model_part.Nodes:
                if control_type == "shape":
                    node.SetSolutionStepValue(KOA.D_LINEAR_D_X, [1,1,1])

            controlle_object_index += 1               
        Logger.PrintInfo("LinearResponseFunction", "Time needed for calculating gradients ",round(timer.time() - startTime,2),"s")  