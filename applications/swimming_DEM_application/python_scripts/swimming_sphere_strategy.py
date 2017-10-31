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
        do_search_neighbours =  self.DEM_parameters["do_search_neighbours"].GetBool()
        print('self.DEM_parameters.IntegrationScheme', self.DEM_parameters["IntegrationScheme"].GetString())
        print('self.DEM_parameters.do_search_neighbours', do_search_neighbours)

        if self.DEM_parameters["IntegrationScheme"].GetString() == 'Verlet_Velocity':
            self.cplusplus_strategy = IterativeSolverStrategy(self.settings, self.max_delta_time, self.n_step_search, self.safety_factor,
                                                              self.delta_option, self.creator_destructor, self.dem_fem_search,
                                                              self.time_integration_scheme, self.search_strategy, do_search_neighbours)

        elif self.DEM_parameters["IntegrationScheme"].GetString() in {'Hybrid_Bashforth', 'TerminalVelocityScheme'}:
            self.cplusplus_strategy = AdamsBashforthStrategy(self.settings, self.max_delta_time, self.n_step_search, self.safety_factor,
                                                              self.delta_option, self.creator_destructor, self.dem_fem_search,
                                                              self.time_integration_scheme, self.search_strategy, do_search_neighbours)

        else:
            self.cplusplus_strategy = ExplicitSolverStrategy(self.settings, self.max_delta_time, self.n_step_search, self.safety_factor,
                                                             self.delta_option, self.creator_destructor, self.dem_fem_search,
                                                             self.time_integration_scheme, self.search_strategy, do_search_neighbours)

    def GetSchemeInstance(self, class_name):
        return globals().get(class_name)()
