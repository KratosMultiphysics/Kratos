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

        phase_input_file_name = parameters["phase_input_file_name"].GetString()
        times_input_file_name = parameters["times_input_file_name"].GetString()
        ambient_input_file_name = parameters["ambient_input_file_name"].GetString()
        
        self.table_phase = PiecewiseLinearTable()
        with open(phase_input_file_name,'r') as file_name1:
            for j, line in enumerate(file_name1):
                file_1 = line.split(" ")
                if (len(file_1)) > 1: 
                    self.table_phase.AddRow(float(file_1[0]), float(file_1[1]))

        self.table_times = PiecewiseLinearTable()
        with open(times_input_file_name,'r') as file_name2:
            for j, line in enumerate(file_name2):
                file_2 = line.split(" ")
                if (len(file_2)) > 1: 
                    self.table_times.AddRow(float(file_2[0]), float(file_2[1]))
        
        self.table_ambient = PiecewiseLinearTable()
        with open(ambient_input_file_name,'r') as file_name3:
            for j, line in enumerate(file_name3):
                file_3 = line.split(" ")
                if (len(file_3)) > 1: 
                    self.table_ambient.AddRow(float(file_3[0]), float(file_3[1]))


        # Construct the utility
        self.Construction = ConstructionUtility(self.mechanical_model_part,self.thermal_model_part, self.table_phase, self.table_times, self.table_ambient, parameters)
    
    def Initialize(self):
        self.Construction.Initialize()

    def InitializeSolutionStep(self):
        self.Construction.InitializeSolutionStep()

    def AfterOutputStep(self):
        self.Construction.AfterOutputStep()