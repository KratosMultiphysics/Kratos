from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
from KratosMultiphysics import *
from KratosMultiphysics.DEMApplication import *
from KratosMultiphysics.SwimmingDEMApplication import *

import sphere_strategy
BaseStrategy = sphere_strategy.ExplicitStrategy

class SwimmingStrategy(BaseStrategy):

    def Initialize(self):
        BaseStrategy.Initialize(self)
        BaseStrategy.SetVariablesAndOptions(self)

        if (self.Parameters.IntegrationScheme == 'Verlet_Velocity'):
            self.cplusplus_strategy = IterativeSolverStrategy(self.settings, self.max_delta_time, self.n_step_search, self.safety_factor,
                                                              self.delta_option, self.creator_destructor, self.dem_fem_search,
                                                              self.time_integration_scheme, self.search_strategy)
            
        elif (self.Parameters.IntegrationScheme == 'Hybrid_Bashforth'):
            self.cplusplus_strategy = AdamsBashforthStrategy(self.settings, self.max_delta_time, self.n_step_search, self.safety_factor,
                                                              self.delta_option, self.creator_destructor, self.dem_fem_search,
                                                              self.time_integration_scheme, self.search_strategy)
        else:
            self.cplusplus_strategy = ExplicitSolverStrategy(self.settings, self.max_delta_time, self.n_step_search, self.safety_factor,
                                                             self.delta_option, self.creator_destructor, self.dem_fem_search,
                                                             self.time_integration_scheme, self.search_strategy)

        self.CheckMomentumConservation()

        self.cplusplus_strategy.Initialize()  # Calls the cplusplus_strategy (C++) Initialize function (initializes all elements and performs other necessary tasks before starting the time loop in Python)
