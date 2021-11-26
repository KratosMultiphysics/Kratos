# This script contains an algorithm that models fluid-particle interaction.
# It combines two parts: a FEM model for the fluid and a DEM model for the particles.
# It has been conceived by adding the DEM part and the interaction on top of an original fluid-only script (see kratos/applications/FluidDynamicswApplication)
# Some parts of the original fluid script have been kept practically untouched and are clearly marked.
# Whenever a minor modification has been made on one of these parts, the corresponding line is indicated with a comment: # MOD.

from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
import warnings

# Kratos
from KratosMultiphysics import *

try:
    from KratosMultiphysics.ExternalSolversApplication import *
except ImportError:
    warnings.warn('Package ExternalSolversApplication could not be loaded. Make sure to compile it if needed.', ImportWarning)

from KratosMultiphysics.DEMApplication import *
from KratosMultiphysics.FluidDynamicsApplication import *
from KratosMultiphysics.SwimmingDEMApplication import *

from KratosMultiphysics.SwimmingDEMApplication.swimming_DEM_analysis import SwimmingDEMAnalysis

class SwimmingDEMAnalysisWithFlush(SwimmingDEMAnalysis):
    def __init__(self, model, algorithm = None, parameters=Parameters("{}")):
        with open('ProjectParameters.json','r') as parameter_file:
            parameters = Parameters(parameter_file.read())
        super(SwimmingDEMAnalysisWithFlush, self).__init__(model, parameters)

    def __enter__ (self):
        return self

    def __exit__(self, exception_type, exception_value, traceback):
        pass

if __name__=="__main__":
    model = Model()
    simulation = SwimmingDEMAnalysisWithFlush(model=model)
    simulation.Run()
