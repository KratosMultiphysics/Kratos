from KratosMultiphysics import *
from KratosMultiphysics.DEMApplication import *

# Thermal DEM Application
import KratosMultiphysics.DEMApplication.sphere_strategy as SolverStrategy

BaseExplicitStrategy = SolverStrategy.ExplicitStrategy

class ExplicitStrategy(BaseExplicitStrategy):

    def __init__(self, all_model_parts, creator_destructor, dem_fem_search, DEM_parameters, procedures):

        BaseExplicitStrategy.__init__(self, all_model_parts, creator_destructor, dem_fem_search, DEM_parameters, procedures)

        # SIMULATION FLAGS

    def AddAdditionalVariables(self, model_part, DEM_parameters):

        BaseExplicitStrategy.AddAdditionalVariables(self, model_part, DEM_parameters)

        model_part.AddNodalSolutionStepVariable(TEMPERATURE)
        model_part.AddNodalSolutionStepVariable(HEATFLUX)



