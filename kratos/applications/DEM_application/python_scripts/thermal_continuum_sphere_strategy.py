from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
from KratosMultiphysics import *
from KratosMultiphysics.DEMApplication import *


# Thermal DEM Application

import continuum_sphere_strategy as SolverStrategy

BaseExplicitStrategy = SolverStrategy.ExplicitStrategy

class ExplicitStrategy(BaseExplicitStrategy):   
   
    def __init__(self, model_part, fem_model_part, cluster_model_part, inlet_model_part, creator_destructor, dem_fem_search, scheme, Param, procedures):

        BaseExplicitStrategy.__init__(self, model_part, fem_model_part, cluster_model_part, inlet_model_part, creator_destructor, dem_fem_search, scheme, Param, procedures)

        # SIMULATION FLAGS  

    def AddAdditionalVariables(self, model_part, Param):
        
        BaseExplicitStrategy.AddAdditionalVariables(self, model_part, Param)

        model_part.AddNodalSolutionStepVariable(TEMPERATURE) 
        model_part.AddNodalSolutionStepVariable(HEATFLUX) 
        
        
        