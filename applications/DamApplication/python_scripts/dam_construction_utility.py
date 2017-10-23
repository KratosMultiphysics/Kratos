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
                self.number_iter = j

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

        time = self.mechanical_model_part.ProcessInfo[TIME]
        delta_time = self.mechanical_model_part.ProcessInfo[DELTA_TIME]
        step = int(time/delta_time)-1

        with open('thermal_parts.txt','r') as file_name4:
            it_t=(thermal_linea for i,thermal_linea in enumerate(file_name4) if i==step)
            for thermal_linea in it_t:
                thermal_name = thermal_linea.rstrip('\n')
        
        with open('mechanical_parts.txt','r') as file_name5:
            it_m=(mechanical_linea for j,mechanical_linea in enumerate(file_name5) if j==step)
            for mechanical_linea in it_m:
                mechanical_name = mechanical_linea.rstrip('\n')

        self.Construction.InitializeSolutionStep(thermal_name,mechanical_name)

        for i in range(self.number_iter+1):
            phase_time = time - self.table_times.GetValue(i)
            if (phase_time>0.0):
                with open('thermal_parts.txt','r') as file_name4:
                    it_t=(thermal_linea for w,thermal_linea in enumerate(file_name4) if w==i)
                    for thermal_linea in it_t:
                        thermal_name = thermal_linea.rstrip('\n')

                self.Construction.ActiveHeatFlux(thermal_name,int(self.table_phase.GetValue(i)),phase_time)

    def AfterOutputStep(self):
        self.Construction.AfterOutputStep()