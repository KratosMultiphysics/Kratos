from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
from KratosMultiphysics import *
from KratosMultiphysics.DEMApplication import *


# Thermal DEM Application

import continuum_sphere_strategy as SolverStrategy

BaseExplicitStrategy = SolverStrategy.ExplicitStrategy

class ExplicitStrategy(BaseExplicitStrategy):   

    def __init__(self, all_model_parts, creator_destructor, dem_fem_search, DEM_parameters, procedures):

        BaseExplicitStrategy.__init__(self, all_model_parts, creator_destructor, dem_fem_search, DEM_parameters, procedures)

        # SIMULATION FLAGS  

    def AddAdditionalVariables(self, spheres_model_part, DEM_parameters):
        
        BaseExplicitStrategy.AddAdditionalVariables(self, spheres_model_part, DEM_parameters)

        spheres_model_part.AddNodalSolutionStepVariable(TEMPERATURE) 
        spheres_model_part.AddNodalSolutionStepVariable(HEATFLUX) 
        
        
        
