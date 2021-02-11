# This script contains an algorithm that models fluid-particle interaction.
# It combines two parts: a FEM model for the fluid and a DEM model for the particles.
# It has been conceived by adding the DEM part and the interaction on top of an original fluid-only script (see kratos/applications/FluidDynamicswApplication)
# Some parts of the original fluid script have been kept practically untouched and are clearly marked.
# Whenever a minor modification has been made on one of these parts, the corresponding line is indicated with a comment: # MOD.

# Kratos
import KratosMultiphysics
from KratosMultiphysics import Model, Parameters
import KratosMultiphysics.ExternalSolversApplication
import KratosMultiphysics.FluidDynamicsApplication
import KratosMultiphysics.DEMApplication
import KratosMultiphysics.SwimmingDEMApplication
from KratosMultiphysics.SwimmingDEMApplication.swimming_DEM_analysis import SwimmingDEMAnalysis

class SwimmingDEMAnalysisWithFlush(SwimmingDEMAnalysis):
    def __init__(self, model, algorithm = None, parameters=Parameters("{}")):
        with open('ProjectParameters.json','r') as parameter_file:
            parameters = Parameters(parameter_file.read())
        super().__init__(model, parameters)

    def __enter__ (self):
        return self

    def __exit__(self, exception_type, exception_value, traceback):
        pass

if __name__=="__main__":
    model = Model()
    simulation = SwimmingDEMAnalysisWithFlush(model=model)
    simulation.Run()
