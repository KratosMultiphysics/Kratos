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

        if name == 'TerminalVelocityScheme':
            class_name = 'TerminalVelocityScheme'

        return class_name

    def CreateCPlusPlusStrategy(self):
        self.SetVariablesAndOptions()
        print('self.Parameters.IntegrationScheme',self.Parameters.IntegrationScheme)
        print('self.Parameters.do_search_neighbours',self.Parameters.do_search_neighbours)

        if self.Parameters.IntegrationScheme == 'Verlet_Velocity':
            self.cplusplus_strategy = IterativeSolverStrategy(self.settings, self.max_delta_time, self.n_step_search, self.safety_factor,
                                                              self.delta_option, self.creator_destructor, self.dem_fem_search,
                                                              self.time_integration_scheme, self.search_strategy, self.Parameters.do_search_neighbours)

        elif self.Parameters.IntegrationScheme in {'Hybrid_Bashforth', 'TerminalVelocityScheme'}:
            self.cplusplus_strategy = AdamsBashforthStrategy(self.settings, self.max_delta_time, self.n_step_search, self.safety_factor,
                                                              self.delta_option, self.creator_destructor, self.dem_fem_search,
                                                              self.time_integration_scheme, self.search_strategy, self.Parameters.do_search_neighbours)

        else:
            self.cplusplus_strategy = ExplicitSolverStrategy(self.settings, self.max_delta_time, self.n_step_search, self.safety_factor,
                                                             self.delta_option, self.creator_destructor, self.dem_fem_search,
                                                             self.time_integration_scheme, self.search_strategy, self.Parameters.do_search_neighbours)


