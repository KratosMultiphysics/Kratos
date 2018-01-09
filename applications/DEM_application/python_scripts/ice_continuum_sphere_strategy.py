from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
from KratosMultiphysics import *
from KratosMultiphysics.DEMApplication import *

# Ice Continuum Strategy

import continuum_sphere_strategy as SolverStrategy

BaseExplicitStrategy = SolverStrategy.ExplicitStrategy

class ExplicitStrategy(BaseExplicitStrategy):   
   
    def __init__(self, all_model_parts, creator_destructor, dem_fem_search, DEM_parameters, procedures):

        BaseExplicitStrategy.__init__(self, all_model_parts, creator_destructor, dem_fem_search, DEM_parameters, procedures)

    def AddAdditionalVariables(self, spheres_model_part, DEM_parameters):
        
        BaseExplicitStrategy.AddAdditionalVariables(self, spheres_model_part, DEM_parameters)

        # Add the necessary variables for the ice_strategy here:
        # spheres_model_part.AddNodalSolutionStepVariable(VARIABLE1)
        # ...
        
        
        
