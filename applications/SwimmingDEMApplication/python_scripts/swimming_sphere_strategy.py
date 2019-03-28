from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
from KratosMultiphysics import *
from KratosMultiphysics.DEMApplication import *
from KratosMultiphysics.SwimmingDEMApplication import *

import sphere_strategy
BaseStrategy = sphere_strategy.ExplicitStrategy

class SwimmingStrategy(BaseStrategy):
    def __init__(self, all_model_parts, creator_destructor, dem_fem_search, parameters, procedures):
        self.project_parameters = parameters
        super(SwimmingStrategy, self).__init__(all_model_parts, creator_destructor, dem_fem_search, parameters['dem_parameters'], procedures)

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
        do_search_neighbours =  self.project_parameters["do_search_neighbours"].GetBool()
        solver_settings = self.DEM_parameters["solver_settings"]

        if self.DEM_parameters["TranslationalIntegrationScheme"].GetString() == 'Verlet_Velocity':
            self.cplusplus_strategy = IterativeSolverStrategy(self.settings, self.max_delta_time, self.n_step_search, self.safety_factor,
                                                              self.delta_option, self.creator_destructor, self.dem_fem_search,
                                                              self.search_strategy, solver_settings, do_search_neighbours)

        elif self.DEM_parameters["TranslationalIntegrationScheme"].GetString() in {'Hybrid_Bashforth', 'TerminalVelocityScheme'}:
            self.cplusplus_strategy = AdamsBashforthStrategy(self.settings, self.max_delta_time, self.n_step_search, self.safety_factor,
                                                              self.delta_option, self.creator_destructor, self.dem_fem_search,
                                                              self.search_strategy, solver_settings, do_search_neighbours)

        else:
            self.cplusplus_strategy = ExplicitSolverStrategy(self.settings, self.max_delta_time, self.n_step_search, self.safety_factor,
                                                             self.delta_option, self.creator_destructor, self.dem_fem_search,
                                                             self.search_strategy, solver_settings, do_search_neighbours)

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

    def GetHydrodynamicLawParametersIfItExists(self, properties):
        if self.project_parameters.Has('properties'):
            for p in self.project_parameters["properties"]:
                if p['properties_id'].GetInt() == int(properties.Id) and p.Has('hydrodynamic_law_parameters'):
                    return p['hydrodynamic_law_parameters']
        return None

    @staticmethod
    def CreateHydrodynamicLaw(properties, hydrodynamic_law_parameters):

        hydrodynamic_name = hydrodynamic_law_parameters['name'].GetString()
        HydrodynamicInteractionLaw = globals().get(hydrodynamic_name)(properties, hydrodynamic_law_parameters)

        if hydrodynamic_law_parameters.Has('buoyancy_parameters'):
            buoyancy_parameters = hydrodynamic_law_parameters['buoyancy_parameters']
            buoyancy_name = buoyancy_parameters['name'].GetString()
            if not buoyancy_name == 'default':
                buoyancy_law = globals().get(buoyancy_name)(buoyancy_parameters)
                HydrodynamicInteractionLaw.SetBuoyancyLaw(buoyancy_law)

        if hydrodynamic_law_parameters.Has('inviscid_force_parameters'):
            inviscid_force_parameters = hydrodynamic_law_parameters['inviscid_force_parameters']
            inviscid_force_name = inviscid_force_parameters['name'].GetString()
            if not inviscid_force_name == 'default':
                inviscid_force_law = globals().get(inviscid_force_name)(inviscid_force_parameters)
                HydrodynamicInteractionLaw.SetInviscidForceLaw(inviscid_force_law)

        if hydrodynamic_law_parameters.Has('drag_parameters'):
            drag_parameters = hydrodynamic_law_parameters['drag_parameters']
            drag_name = drag_parameters['name'].GetString()
            if not drag_name == 'default':
                drag_law = globals().get(drag_name)(drag_parameters)
                HydrodynamicInteractionLaw.SetDragLaw(drag_law)

        if hydrodynamic_law_parameters.Has('history_force_parameters'):
            history_force_parameters = hydrodynamic_law_parameters['history_force_parameters']
            history_force_name = history_force_parameters['name'].GetString()
            if not history_force_name == 'default':
                history_force_law = globals().get(history_force_name)(history_force_parameters)
                HydrodynamicInteractionLaw.SetHistoryForceLaw(history_force_law)

        if hydrodynamic_law_parameters.Has('vorticity_induced_lift_parameters'):
            vorticity_induced_lift_parameters = hydrodynamic_law_parameters['vorticity_induced_lift_parameters']
            vorticity_induced_lift_name = vorticity_induced_lift_parameters['name'].GetString()
            if not vorticity_induced_lift_name == 'default':
                vorticity_induced_lift_law = globals().get(vorticity_induced_lift_name)(vorticity_induced_lift_parameters)
                HydrodynamicInteractionLaw.SetVorticityInducedLiftLaw(vorticity_induced_lift_law)

        if hydrodynamic_law_parameters.Has('rotation_induced_lift_parameters'):
            rotation_induced_lift_parameters = hydrodynamic_law_parameters['rotation_induced_lift_parameters']
            rotation_induced_lift_name = rotation_induced_lift_parameters['name'].GetString()
            if not rotation_induced_lift_name == 'default':
                rotation_induced_lift_law = globals().get(rotation_induced_lift_name)(rotation_induced_lift_parameters)
                HydrodynamicInteractionLaw.SetRotationInducedLiftLaw(rotation_induced_lift_law)

        if hydrodynamic_law_parameters.Has('steady_viscous_torque_parameters'):
            steady_viscous_torque_parameters = hydrodynamic_law_parameters['steady_viscous_torque_parameters']
            steady_viscous_torque_name = steady_viscous_torque_parameters['name'].GetString()
            if not steady_viscous_torque_name == 'default':
                steady_viscous_torque_law = globals().get(steady_viscous_torque_name)(steady_viscous_torque_parameters)
                HydrodynamicInteractionLaw.SetSteadyViscousTorqueLaw(steady_viscous_torque_law)

        HydrodynamicInteractionLaw.SetHydrodynamicInteractionLawInProperties(properties, True)

    def ModifyProperties(self, properties, param = 0):

        super(SwimmingStrategy,self).ModifyProperties(properties, param)

        hydrodynamic_law_parameters = self.GetHydrodynamicLawParametersIfItExists(properties)
        if hydrodynamic_law_parameters:
            SwimmingStrategy.CreateHydrodynamicLaw(properties, hydrodynamic_law_parameters)

        if not param:
            if not properties.Has(PARTICLE_SPHERICITY):
                properties[PARTICLE_SPHERICITY] = 1.0