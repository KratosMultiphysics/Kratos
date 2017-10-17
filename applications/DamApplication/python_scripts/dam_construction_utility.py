from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
#importing the Kratos Library
from KratosMultiphysics import *
from KratosMultiphysics.SolidMechanicsApplication import *
from KratosMultiphysics.ExternalSolversApplication import *
from KratosMultiphysics.DamApplication import *

class DamConstructionUtility:

    def __init__(self,mechanical_model_part, thermal_model_part, parameters):

        self.mechanical_model_part = mechanical_model_part
        self.thermal_model_part = thermal_model_part

        # Construct the utility
        self.Construction = ConstructionUtility(self.mechanical_model_part,self.thermal_model_part, parameters)
    
    def Initialize(self):
        self.Construction.Initialize()

    def InitializeSolutionStep(self):
        self.Construction.InitializeSolutionStep()

    def AfterOutputStep(self):
        self.Construction.AfterOutputStep()