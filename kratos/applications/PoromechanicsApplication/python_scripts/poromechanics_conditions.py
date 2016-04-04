from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# importing the Kratos Library
from KratosMultiphysics import *
from KratosMultiphysics.PoromechanicsApplication import *
CheckForPreviousImport()
from math import *


class PoromechanicsConditions:

    def __init__(self, delta_time, ConditionsOptions):

        self.delta_time = delta_time

        self.imposed_displacement = ConditionsOptions.Imposed_Displacement
        self.imposed_pressure = ConditionsOptions.Imposed_Pressure
        self.imposed_force = ConditionsOptions.Imposed_Force
        self.imposed_faceload = ConditionsOptions.Imposed_FaceLoad
        self.imposed_normalload = ConditionsOptions.Imposed_NormalLoad
        self.imposed_tangentialload = ConditionsOptions.Imposed_TangentialLoad
        self.imposed_normalfluidflux = ConditionsOptions.Imposed_NormalFluidFlux
        
        self.listofprocesses = []


    def Initialize(self, model_part):
        
        if(self.imposed_displacement=="Table_Interpolation"):
            self.listofprocesses.append(DisplacementTableInterpolationProcess(model_part))
        
        if(self.imposed_pressure=="Table_Interpolation"):
            self.listofprocesses.append(PressureTableInterpolationProcess(model_part))

        if(self.imposed_force=="Table_Interpolation"):
            self.listofprocesses.append(ForceTableInterpolationProcess(model_part))
            
        if(self.imposed_faceload=="Table_Interpolation"):
            self.listofprocesses.append(FaceLoadTableInterpolationProcess(model_part))
                        
        if(self.imposed_normalload=="Table_Interpolation"):
            self.listofprocesses.append(NormalLoadTableInterpolationProcess(model_part))

        if(self.imposed_tangentialload=="Table_Interpolation"):
            self.listofprocesses.append(TangentialLoadTableInterpolationProcess(model_part))
            
        if(self.imposed_normalfluidflux=="Table_Interpolation"):
            self.listofprocesses.append(NormalFluxTableInterpolationProcess(model_part))


    def UpdateImposedConditions(self, model_part, current_step):
        
        for process in self.listofprocesses:
            process.Execute()
