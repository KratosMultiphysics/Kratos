from KratosMultiphysics import *
from KratosMultiphysics.DEMApplication import *

# Ice Continuum Strategy

import KratosMultiphysics.DEMApplication.continuum_sphere_strategy as SolverStrategy

BaseExplicitStrategy = SolverStrategy.ExplicitStrategy

class ExplicitStrategy(BaseExplicitStrategy):

    def __init__(self, all_model_parts, creator_destructor, dem_fem_search, DEM_parameters, procedures):

        BaseExplicitStrategy.__init__(self, all_model_parts, creator_destructor, dem_fem_search, DEM_parameters, procedures)

    def AddAdditionalVariables(self, spheres_model_part, DEM_parameters):

        BaseExplicitStrategy.AddAdditionalVariables(self, spheres_model_part, DEM_parameters)

        # Add the necessary variables for the ice_strategy here:
        # spheres_model_part.AddNodalSolutionStepVariable(VARIABLE1)
        # ...



