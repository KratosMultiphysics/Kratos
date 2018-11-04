from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
from KratosMultiphysics import *
from KratosMultiphysics.DEMApplication import *
from KratosMultiphysics.SwimmingDEMApplication import *

import sphere_strategy
BaseStrategy = sphere_strategy.ExplicitStrategy

class SwimmingStrategy(BaseStrategy):

    def TranslationalIntegrationSchemeTranslator(self, name):
        class_name = BaseStrategy.TranslationalIntegrationSchemeTranslator(self, name)

        if name == 'Hybrid_Bashforth':
            class_name = 'HybridBashforthScheme'
        if name == 'TerminalVelocityScheme':
            class_name = 'TerminalVelocityScheme'

        return class_name

    def RotationalIntegrationSchemeTranslator(self, name_translational, name_rotational):
        class_name = BaseStrategy.RotationalIntegrationSchemeTranslator(self, name_translational, name_rotational)

        if name_rotational == 'Direct_Integration':
            if name_translational == 'Hybrid_Bashforth':
                class_name = 'HybridBashforthScheme'
            elif name_translational == 'TerminalVelocityScheme':
                class_name = 'TerminalVelocityScheme'
        elif name_rotational == 'Runge_Kutta':
            class_name = 'RungeKuttaScheme'
        elif name_rotational == 'Quaternion_Integration':
            class_name = 'QuaternionIntegrationScheme'

        return class_name

    def CreateCPlusPlusStrategy(self):
        self.SetVariablesAndOptions()
        do_search_neighbours =  self.DEM_parameters["do_search_neighbours"].GetBool()
        strategy_parameters = self.DEM_parameters["strategy_parameters"]

        if self.DEM_parameters["TranslationalIntegrationScheme"].GetString() == 'Verlet_Velocity':
            self.cplusplus_strategy = IterativeSolverStrategy(self.settings, self.max_delta_time, self.n_step_search, self.safety_factor,
                                                              self.delta_option, self.creator_destructor, self.dem_fem_search,
                                                              self.search_strategy, strategy_parameters, do_search_neighbours)

        elif self.DEM_parameters["TranslationalIntegrationScheme"].GetString() in {'Hybrid_Bashforth', 'TerminalVelocityScheme'}:
            self.cplusplus_strategy = AdamsBashforthStrategy(self.settings, self.max_delta_time, self.n_step_search, self.safety_factor,
                                                              self.delta_option, self.creator_destructor, self.dem_fem_search,
                                                              self.search_strategy, strategy_parameters, do_search_neighbours)

        else:
            self.cplusplus_strategy = ExplicitSolverStrategy(self.settings, self.max_delta_time, self.n_step_search, self.safety_factor,
                                                             self.delta_option, self.creator_destructor, self.dem_fem_search,
                                                             self.search_strategy, strategy_parameters, do_search_neighbours)

    def GetTranslationalSchemeInstance(self, class_name):
         if not class_name == 'NewmarkBetaScheme':
             return globals().get(class_name)()
         else:
             return globals().get(class_name)(0.5,0.25)

    def GetRotationalSchemeInstance(self, class_name):
         if not class_name == 'NewmarkBetaScheme':
             return globals().get(class_name)()
         else:
             return globals().get(class_name)(0.5,0.25)

    def ModifyProperties(self, properties, param = 0):

        super(SwimmingStrategy,self).ModifyProperties(properties, param)

        if not param:
            if not properties.Has(PARTICLE_SPHERICITY):
                properties[PARTICLE_SPHERICITY] = 1.0
