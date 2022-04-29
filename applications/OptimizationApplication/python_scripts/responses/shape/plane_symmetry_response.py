import KratosMultiphysics as KM
from KratosMultiphysics import Parameters, Logger
import KratosMultiphysics.OptimizationApplication as KOA
from KratosMultiphysics.OptimizationApplication.responses.base_response import BaseResponseFunction
import KratosMultiphysics.StructuralMechanicsApplication as KSM

import time as timer
import numpy as np

class PlaneSymmetryResponseFunction(BaseResponseFunction):

    def __init__(self,response_name, response_settings,model):

        self.type = "plane_symmetry"
        self.variable = "PLANE_SYMMETRY"
        super().__init__(response_name, response_settings, model)        
    
        self.supported_control_types = ["shape"]
        self.gradients_variables = {"shape":"D_PLANE_SYMMETRY_D_X"}

        if len(self.evaluated_model_parts) != len(self.controlled_model_parts) or set(self.evaluated_model_parts) != set(self.controlled_model_parts):
            raise RuntimeError("PlaneSymmetryResponseFunction: 'evaluated_objects' of response '{}' must be the same as controlled_objects !".format(self.name)) 

        for control_type in self.control_types:
            if not control_type in self.supported_control_types:
                raise RuntimeError("PlaneSymmetryResponseFunction: type {} in 'control_types' of response '{}' is not supported, supported types are {}  !".format(control_type,self.name,self.supported_control_types)) 

        
        root_model_part_name = self.evaluated_model_parts[0].split(".")[0]
        for evaluated_model_part in self.evaluated_model_parts:
            if evaluated_model_part.split(".")[0] != root_model_part_name:
                raise RuntimeError("PlaneSymmetryResponseFunction: evaluated_model_parts of plane symmetry response must have the same root model part !")

        self.root_model_part = self.model.GetModelPart(root_model_part_name)
        # add vars
        for control_type in self.control_types:
            if control_type == "shape":
                self.root_model_part.AddNodalSolutionStepVariable(KM.KratosGlobals.GetVariable(self.gradients_variables[control_type]))
                self.root_model_part.AddNodalSolutionStepVariable(KOA.NEAREST_NEIGHBOUR_POINT)
                self.root_model_part.AddNodalSolutionStepVariable(KOA.NEAREST_NEIGHBOUR_DIST)

        self.plane_symmetry_response = KOA.PlaneSymmetry(response_name,model,self.response_settings)

    def GetVariableName(self):
        return  self.variable

    def GetGradientsVariablesName(self):
        return self.gradients_variables

    def GetGradientVariableNameForType(self,control_type, raise_error=True):
        if raise_error:
            if not control_type in self.supported_control_types:
                raise RuntimeError("PlaneSymmetryResponseFunction: type {} in 'control_types' of response '{}' is not supported, supported types are {}  !".format(control_type,self.name,self.supported_control_types)) 

        return self.gradients_variables[control_type]

    def Initialize(self):
        super().Initialize()
        # self.plane_symmetry_response.Initialize()

    def CalculateValue(self):
        Logger.PrintInfo("PlaneSymmetryResponseFunction:CalculateValue: Starting value calculation for response ", self.name)
        startTime = timer.time()
        self.value = 2.0
        # self.plane_symmetry_response.CalculateValue()
        Logger.PrintInfo("PlaneSymmetryResponseFunction:CalculateValue: Time needed for calculating value ",round(timer.time() - startTime,2),"s")        
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
                    raise RuntimeError("PlaneSymmetryResponseFunction:CalculateGradientsForTypesAndObjects: control type {} of control object {} is not in the control_types of response {}".format(control_types[itr],controlled_object,self.name))

        Logger.PrintInfo("PlaneSymmetryResponseFunction", "Starting ", control_types," gradients calculation of response ", self.name," for ",controlled_objects)
        startTime = timer.time()

        # self.plane_symmetry_response.CalculateGradient()
               
        Logger.PrintInfo("PlaneSymmetryResponseFunction", "Time needed for calculating gradients ",round(timer.time() - startTime,2),"s")  