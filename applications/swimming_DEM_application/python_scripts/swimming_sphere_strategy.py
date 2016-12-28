from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
from KratosMultiphysics import *
from KratosMultiphysics.DEMApplication import *
from KratosMultiphysics.SwimmingDEMApplication import *

import sphere_strategy
BaseStrategy = sphere_strategy.ExplicitStrategy

class SwimmingStrategy(BaseStrategy):

    def IntegrationSchemeTranslator(self, name):
        class_name = BaseStrategy.IntegrationSchemeTranslator(self, name)

        if name == 'Hybrid_Bashforth':
            class_name = 'HybridBashforthScheme'

        return class_name

    def GetSchemeInstance(self, class_name): # parent counterpart must not be called due to different 'globals()'
        return globals().get(class_name)()

    def Initialize(self):
        BaseStrategy.Initialize(self)
        BaseStrategy.SetVariablesAndOptions(self)

        self.CheckMomentumConservation()

        self.cplusplus_strategy.Initialize()  # Calls the cplusplus_strategy (C++) Initialize function (initializes all elements and performs other necessary tasks before starting the time loop in Python)

    def CreateCPlusPlusStrategy(self):
        self.SetVariablesAndOptions()
        print('self.Parameters.IntegrationScheme',self.Parameters.IntegrationScheme)
        print('self.Parameters.do_search_neighbours',self.Parameters.do_search_neighbours)

        if (self.Parameters.IntegrationScheme == 'Verlet_Velocity'):
            self.cplusplus_strategy = IterativeSolverStrategy(self.settings, self.max_delta_time, self.n_step_search, self.safety_factor,
                                                              self.delta_option, self.creator_destructor, self.dem_fem_search,
                                                              self.time_integration_scheme, self.search_strategy, self.Parameters.do_search_neighbours)

        elif (self.Parameters.IntegrationScheme == 'Hybrid_Bashforth'):
            self.cplusplus_strategy = AdamsBashforthStrategy(self.settings, self.max_delta_time, self.n_step_search, self.safety_factor,
                                                              self.delta_option, self.creator_destructor, self.dem_fem_search,
                                                              self.time_integration_scheme, self.search_strategy, self.Parameters.do_search_neighbours)

        else:
            self.cplusplus_strategy = ExplicitSolverStrategy(self.settings, self.max_delta_time, self.n_step_search, self.safety_factor,
                                                             self.delta_option, self.creator_destructor, self.dem_fem_search,
                                                             self.time_integration_scheme, self.search_strategy, self.Parameters.do_search_neighbours)


