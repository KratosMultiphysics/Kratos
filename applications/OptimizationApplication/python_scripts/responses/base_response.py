# importing the Kratos Library
from numpy import gradient
import KratosMultiphysics as KM
from KratosMultiphysics import Parameters, Logger
import KratosMultiphysics.OptimizationApplication as KOA

import time as timer
import numpy as np

# ==============================================================================
class BaseResponseFunction:

    def __init__(self,response_name, response_settings, model, response_analysis=None):
       
        self.name = response_name
        self.response_settings = response_settings
        self.model = model
        self.analysis = response_analysis
        self.analysis_model_part = None
        if not response_analysis == None:
            self.analysis_model_part = self.analysis._GetSolver().GetComputingModelPart()

        self.evaluated_model_parts = response_settings["evaluated_objects"].GetStringArray()
        self.controlled_model_parts = response_settings["controlled_objects"].GetStringArray()
        self.control_types = response_settings["control_types"].GetStringArray()  

    def Initialize(self):

        for evaluated_model_part in self.evaluated_model_parts:
            if not self.analysis == None:
                evaluated_model_part_splitted = evaluated_model_part.split(".")
                if not evaluated_model_part_splitted[0] == self.analysis_model_part.Name:
                    raise RuntimeError("BaseResponseFunction:Initialize: root evaluated_model_part {} of response '{}' is not the analysis model!".format(evaluated_model_part_splitted[0],self.name))
            if not self.model.HasModelPart(evaluated_model_part): 
                raise RuntimeError("BaseResponseFunction:Initialize: evaluated_model_part {} of response '{}' does not exist!".format(evaluated_model_part,self.name))

        for controlled_model_part in self.controlled_model_parts:
            if not self.model.HasModelPart(controlled_model_part): 
                raise RuntimeError("BaseResponseFunction:Initialize: controlled_model_part {} of response '{}' does not exist!".format(controlled_model_part,self.name))

    def CalculateValue(self):
        raise RuntimeError("BaseResponseFunction:CalculateValue: Not implemeted ! ")

    def GetValue(self):
        raise RuntimeError("BaseResponseFunction:GetValue: Not implemeted ! ")

    def CalculateGradients(self):
        raise RuntimeError("BaseResponseFunction:CalculateGradients: Not implemeted ! ")

    def CalculateGradientsForTypeAndObjects(self,control_type,controlled_objects,raise_error=True):
        raise RuntimeError("BaseResponseFunction:CalculateGradientsForTypeAndObjects: Not implemeted ! ") 

    def GetGradients(self):
        raise RuntimeError("BaseResponseFunction:GetGradients: Not implemeted ! ")

    def GetType(self):
        raise RuntimeError("BaseResponseFunction:GetValue: Not implemeted ! ")        

    def GetVariableName(self):
        raise RuntimeError("BaseResponseFunction:GetVariableName: Not implemeted ! ") 

    def GetGradientsVariablesName(self):
        raise RuntimeError("BaseResponseFunction:GetGradientsVariablesName: Not implemeted ! ")

    def GetGradientVariableNameForType(self,control_type):
        raise RuntimeError("BaseResponseFunction:GetGradientVariableNameForType: Not implemeted ! ") 