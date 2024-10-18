# By Rafael Rangel (2022)
#   This class was brought from the DEMApp and
#   is intended to be used with the sintering particle.

# Imports
from   KratosMultiphysics import *
from   KratosMultiphysics.DEMApplication import *
from   KratosMultiphysics.ThermalDEMApplication import *
import KratosMultiphysics.DEMApplication.continuum_sphere_strategy as SolverStrategy

# Set base class
BaseStrategy = SolverStrategy.ExplicitStrategy

# Strategy class
class ExplicitStrategy(BaseStrategy):
    #----------------------------------------------------------------------------------------------
    def __init__(self, all_model_parts, creator_destructor, dem_fem_search, DEM_parameters, procedures):
        # Initialize base class
        BaseStrategy.__init__(self, all_model_parts, creator_destructor, dem_fem_search, DEM_parameters, procedures)

    #----------------------------------------------------------------------------------------------
    def AddAdditionalVariables(self, spheres_model_part, DEM_parameters):
        BaseStrategy.AddAdditionalVariables(self, spheres_model_part, DEM_parameters)

        spheres_model_part.AddNodalSolutionStepVariable(TEMPERATURE)
        spheres_model_part.AddNodalSolutionStepVariable(HEATFLUX)
