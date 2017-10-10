from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
#importing the Kratos Library
from KratosMultiphysics import *
from KratosMultiphysics.SolidMechanicsApplication import *
from KratosMultiphysics.ExternalSolversApplication import *
from KratosMultiphysics.DamApplication import *

class DamConstructionUtility:

    def __init__(self,mechanical_model_part, thermal_model_part):

        self.mechanical_model_part = mechanical_model_part
        self.thermal_model_part = thermal_model_part

        self.parameters = Parameters("{}")
        self.parameters.AddEmptyValue("mesh_id").SetInt(0)
        self.parameters.AddEmptyValue("Gravity_Direction").SetString("Z")
        self.parameters.AddEmptyValue("Reservoir_Bottom_Coordinate_in_Gravity_Direction").SetDouble(0.0)
        self.parameters.AddEmptyValue("Height_Dam").SetDouble(1.0)
        self.parameters.AddEmptyValue("Number_of_phases").SetInt(3)

        # Construct the utility
        self.Construction = ConstructionUtility(self.mechanical_model_part,self.thermal_model_part,self.parameters)
    
    def Initialize(self):
        self.Construction.Initialize()

    def InitializeSolutionStep(self):
        self.Construction.InitializeSolutionStep()

    def AfterOutputStep(self):
        self.Construction.AfterOutputStep()