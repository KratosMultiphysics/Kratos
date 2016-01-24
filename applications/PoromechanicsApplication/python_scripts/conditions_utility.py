from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# importing the Kratos Library
from KratosMultiphysics import *
from KratosMultiphysics.PoromechanicsApplication import *
CheckForPreviousImport()
from math import *


class ConditionsUtility:

    def __init__(self, delta_time, ConditionsOptions):

        self.delta_time = delta_time

        self.imposed_displacement = ConditionsOptions.Imposed_Displacement
        self.imposed_pressure = ConditionsOptions.Imposed_Pressure
        self.imposed_pointload = ConditionsOptions.Imposed_PointLoad
        self.imposed_lineload = ConditionsOptions.Imposed_LineLoad
        self.imposed_surfaceload = ConditionsOptions.Imposed_SurfaceLoad
        self.imposed_normalload = ConditionsOptions.Imposed_NormalLoad
        self.imposed_tangentialload = ConditionsOptions.Imposed_TangentialLoad
        self.imposed_normalfluidflux = ConditionsOptions.Imposed_NormalFluidFlux


    def Initialize(self, model_part):

        NonLinearImposedConditions = False

        if(self.imposed_displacement=="Linearly_Incremented"):
            for node in model_part.Nodes:
                ImposedDisplacement = node.GetSolutionStepValue(IMPOSED_DISPLACEMENT)
                Displacement = node.GetSolutionStepValue(DISPLACEMENT)
                if(node.IsFixed(DISPLACEMENT_X)):
                    ImposedDisplacement[0] = Displacement[0]
                    Displacement[0] = 0
                if(node.IsFixed(DISPLACEMENT_Y)):
                    ImposedDisplacement[1] = Displacement[1]
                    Displacement[1] = 0
                if(node.IsFixed(DISPLACEMENT_Z)):
                    ImposedDisplacement[2] = Displacement[2]
                    Displacement[2] = 0
                node.SetSolutionStepValue(IMPOSED_DISPLACEMENT, ImposedDisplacement)
                node.SetSolutionStepValue(DISPLACEMENT, Displacement)
        elif(self.imposed_displacement=="Non-linearly_Modified"):
            NonLinearImposedConditions = True

        if(self.imposed_pressure=="Linearly_Incremented"):
            for node in model_part.Nodes:
                if(node.IsFixed(WATER_PRESSURE)):
                    FluidPressure = node.GetSolutionStepValue(WATER_PRESSURE)
                    ImposedFluidPressure = FluidPressure
                    FluidPressure = 0
                    node.SetSolutionStepValue(IMPOSED_FLUID_PRESSURE, ImposedFluidPressure)
                    node.SetSolutionStepValue(WATER_PRESSURE, FluidPressure)
        elif(self.imposed_pressure=="Non-linearly_Modified"):
            NonLinearImposedConditions = True

        if(self.imposed_pointload=="Linearly_Incremented"):
            for node in model_part.Nodes:
                ImposedPointLoad = node.GetSolutionStepValue(IMPOSED_POINT_LOAD)
                PointLoad = node.GetSolutionStepValue(POINT_LOAD)
                if(node.IsFixed(POINT_LOAD_X)):
                    ImposedPointLoad[0] = PointLoad[0]
                    PointLoad[0] = 0
                if(node.IsFixed(POINT_LOAD_Y)):
                    ImposedPointLoad[1] = PointLoad[1]
                    PointLoad[1] = 0
                if(node.IsFixed(POINT_LOAD_Z)):
                    ImposedPointLoad[2] = PointLoad[2]
                    PointLoad[2] = 0
                node.SetSolutionStepValue(IMPOSED_POINT_LOAD, ImposedPointLoad)
                node.SetSolutionStepValue(POINT_LOAD, PointLoad)
        elif(self.imposed_pointload=="Non-linearly_Modified"):
            NonLinearImposedConditions = True

        if(self.imposed_lineload=="Linearly_Incremented"):
            for node in model_part.Nodes:
                ImposedLineLoad = node.GetSolutionStepValue(IMPOSED_LINE_LOAD)
                LineLoad = node.GetSolutionStepValue(LINE_LOAD)
                if(node.IsFixed(LINE_LOAD_X)):
                    ImposedLineLoad[0] = LineLoad[0]
                    LineLoad[0] = 0
                if(node.IsFixed(LINE_LOAD_Y)):
                    ImposedLineLoad[1] = LineLoad[1]
                    LineLoad[1] = 0
                if(node.IsFixed(LINE_LOAD_Z)):
                    ImposedLineLoad[2] = LineLoad[2]
                    LineLoad[2] = 0
                node.SetSolutionStepValue(IMPOSED_LINE_LOAD, ImposedLineLoad)
                node.SetSolutionStepValue(LINE_LOAD, LineLoad)
        elif(self.imposed_lineload=="Non-linearly_Modified"):
            NonLinearImposedConditions = True

        if(self.imposed_surfaceload=="Linearly_Incremented"):
            for node in model_part.Nodes:
                ImposedSurfaceLoad = node.GetSolutionStepValue(IMPOSED_SURFACE_LOAD)
                SurfaceLoad = node.GetSolutionStepValue(SURFACE_LOAD)
                if(node.IsFixed(SURFACE_LOAD_X)):
                    ImposedSurfaceLoad[0] = SurfaceLoad[0]
                    SurfaceLoad[0] = 0
                if(node.IsFixed(SURFACE_LOAD_Y)):
                    ImposedSurfaceLoad[1] = SurfaceLoad[1]
                    SurfaceLoad[1] = 0
                if(node.IsFixed(SURFACE_LOAD_Z)):
                    ImposedSurfaceLoad[2] = SurfaceLoad[2]
                    SurfaceLoad[2] = 0
                node.SetSolutionStepValue(IMPOSED_SURFACE_LOAD, ImposedSurfaceLoad)
                node.SetSolutionStepValue(SURFACE_LOAD, SurfaceLoad)
        elif(self.imposed_surfaceload=="Non-linearly_Modified"):
            NonLinearImposedConditions = True

        if(self.imposed_normalload=="Linearly_Incremented"):
            for node in model_part.Nodes:
                if(node.IsFixed(NORMAL_CONTACT_STRESS)):
                    NormalStress = node.GetSolutionStepValue(NORMAL_CONTACT_STRESS)
                    ImposedNormalStress = NormalStress
                    NormalStress = 0
                    node.SetSolutionStepValue(IMPOSED_NORMAL_STRESS, ImposedNormalStress)
                    node.SetSolutionStepValue(NORMAL_CONTACT_STRESS, NormalStress)
        elif(self.imposed_normalload=="Non-linearly_Modified"):
            NonLinearImposedConditions = True

        if(self.imposed_tangentialload=="Linearly_Incremented"):
            for node in model_part.Nodes:
                if(node.IsFixed(TANGENTIAL_CONTACT_STRESS)):
                    TangentialStress = node.GetSolutionStepValue(TANGENTIAL_CONTACT_STRESS)
                    ImposedTangentialStress = TangentialStress
                    TangentialStress = 0
                    node.SetSolutionStepValue(IMPOSED_TANGENTIAL_STRESS, ImposedTangentialStress)
                    node.SetSolutionStepValue(TANGENTIAL_CONTACT_STRESS, TangentialStress)
        elif(self.imposed_tangentialload=="Non-linearly_Modified"):
            NonLinearImposedConditions = True

        if(self.imposed_normalfluidflux=="Linearly_Incremented"):
            for node in model_part.Nodes:
                if(node.IsFixed(NORMAL_FLUID_FLUX)):
                    FluidFlux = node.GetSolutionStepValue(NORMAL_FLUID_FLUX)
                    ImposedNormalFluidFlux = FluidFlux
                    FluidFlux = 0
                    node.SetSolutionStepValue(IMPOSED_NORMAL_FLUID_FLUX, ImposedNormalFluidFlux)
                    node.SetSolutionStepValue(NORMAL_FLUID_FLUX, FluidFlux)
        elif(self.imposed_normalfluidflux=="Non-linearly_Modified"):
            NonLinearImposedConditions = True

        return NonLinearImposedConditions


    def UpdateImposedConditions(self, model_part, current_step):
        
        print("Customize your non-linear conditions")
        
        #Soil column problem step load (2D)
        '''
        if(current_step==1):
            for node in model_part.Nodes:
                ImposedLineLoad = node.GetSolutionStepValue(IMPOSED_LINE_LOAD)
                LineLoad = node.GetSolutionStepValue(LINE_LOAD)
                if(node.IsFixed(LINE_LOAD_Y)):
                    ImposedLineLoad[1] = -200
                    LineLoad[1] = 0
                node.SetSolutionStepValue(IMPOSED_LINE_LOAD, ImposedLineLoad)
                node.SetSolutionStepValue(LINE_LOAD, LineLoad)

        if(current_step==6):
            for node in model_part.Nodes:
                ImposedLineLoad = node.GetSolutionStepValue(IMPOSED_LINE_LOAD)
                if(node.IsFixed(LINE_LOAD_Y)):
                    ImposedLineLoad[1] = 0
                node.SetSolutionStepValue(IMPOSED_LINE_LOAD, ImposedLineLoad)
        '''

        #Soil column problem cyclic load (2D)
        '''
        current_time = current_step*self.delta_time
        ImposedLoad = -1e3*(1+0.5*sin(2*current_time))
        for node in model_part.Nodes:
            LineLoad = node.GetSolutionStepValue(LINE_LOAD)
            if(node.IsFixed(LINE_LOAD_Y)):
                LineLoad[1] = ImposedLoad
            node.SetSolutionStepValue(LINE_LOAD, LineLoad)
        '''
