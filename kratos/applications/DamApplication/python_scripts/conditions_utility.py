from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# importing the Kratos Library
from KratosMultiphysics import *
from KratosMultiphysics.SolidMechanicsApplication import *
from KratosMultiphysics.DamApplication import *
CheckForPreviousImport()
from math import *

class ConditionsUtility:


    def __init__(self, delta_time, ConditionsOptions, model_part, time_unit_converter,evolution_type):

        self.time_unit_converter = time_unit_converter
        self.delta_time = delta_time
        self.evolution_type = evolution_type

        self.imposed_displacement = ConditionsOptions.Imposed_Displacement
        self.imposed_temperature = ConditionsOptions.Imposed_Temperature
        self.imposed_pointload = ConditionsOptions.Imposed_PointLoad
        self.imposed_lineload = ConditionsOptions.Imposed_LineLoad
        self.imposed_surfaceload = ConditionsOptions.Imposed_SurfaceLoad
        self.imposed_normalload= ConditionsOptions.Imposed_NormalLoad
        self.imposed_waterload = ConditionsOptions.Imposed_WaterLoad
        self.imposed_Bofang_temperature = ConditionsOptions.Imposed_Bofang_Temperature

        self.listofprocesses=[]
        
    def Initialize(self, model_part):
      
        if(self.imposed_displacement=="Table_Interpolation"):
            self.listofprocesses.append(DisplaTableInterpolationProcess(model_part, self.time_unit_converter))
        
        if(self.imposed_temperature=="Table_Interpolation"):
            self.listofprocesses.append(TemperatureTableInterpolationProcess(model_part, self.time_unit_converter))
            
        if(self.imposed_pointload=="Table_Interpolation"):
            self.listofprocesses.append(PointLoadTableInterpolationProcess(model_part, self.time_unit_converter))
            
        if(self.imposed_lineload=="Table_Interpolation"):
            self.listofprocesses.append(LineLoadTableInterpolationProcess(model_part, self.time_unit_converter))
            
        if(self.imposed_surfaceload=="Table_Interpolation"):
            self.listofprocesses.append(SurfaceLoadTableInterpolationProcess(model_part, self.time_unit_converter))
            
        if(self.imposed_normalload=="Table_Interpolation"):
            self.listofprocesses.append(NormalLoadTableInterpolationProcess(model_part, self.time_unit_converter))
        
        if(self.imposed_waterload=="Table_Evolution_Data"):
            if(self.evolution_type== "Exact"):
                self.listofprocesses.append(ExactWaterEvolutionConditionsLoadProcess(model_part, self.time_unit_converter))
            else:
                self.listofprocesses.append(InterpolationWaterEvolutionConditionsLoadProcess(model_part, self.time_unit_converter))
                
        if(self.imposed_Bofang_temperature=="Table_Evolution_Data"):
            if(self.evolution_type== "Exact"):
                self.listofprocesses.append(ExactBofangEvolutionConditionsTemperatureProcess(model_part, self.time_unit_converter))
            else: 
                self.listofprocesses.append(InterpolationBofangEvolutionConditionsTemperatureProcess(model_part, self.time_unit_converter))

    def UpdateImposedConditions(self, model_part, current_step):
            
        for process in self.listofprocesses:
            process.Execute()
       
