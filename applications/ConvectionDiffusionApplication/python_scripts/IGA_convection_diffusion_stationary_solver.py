import KratosMultiphysics
import numpy as np

# Import base class file
from KratosMultiphysics.ConvectionDiffusionApplication import convection_diffusion_stationary_solver

def CreateSolver(main_model_part, custom_settings):
    return IGAConvectionDiffusionStationarySolver(main_model_part, custom_settings)


class IGAConvectionDiffusionStationarySolver(convection_diffusion_stationary_solver.ConvectionDiffusionStationarySolver):
    print('IGA ci siamo')

    def __init__(self, main_model_part, custom_settings):
        super().__init__(main_model_part, custom_settings)
        KratosMultiphysics.Logger.PrintInfo(self.__class__.__name__, "Construction finished")

    def AddVariables(self):
        super().AddVariables()

    def Initialize(self):
        super().Initialize()

        
    def InitializeSolutionStep(self):
        super().InitializeSolutionStep()
   
    
    def FinalizeSolutionStep(self):
        super().FinalizeSolutionStep()

        

    def Finalize(self):
        super().Finalize()





