import KratosMultiphysics as Kratos
import KratosMultiphysics.DEMApplication as DEM
import KratosMultiphysics.SwimmingDEMApplication as SDEM

from KratosMultiphysics.DEMApplication.sphere_strategy import ExplicitStrategy
BaseStrategy = ExplicitStrategy

class SwimmingStrategy(BaseStrategy):

    @staticmethod
    def SDEMEvaluateString(name):
        return getattr(SDEM, name)

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

        if self.DEM_parameters["TranslationalIntegrationScheme"].GetString() == 'Verlet_Velocity':
            self.cplusplus_strategy = DEM.IterativeSolverStrategy(self.settings, self.max_delta_time, self.n_step_search, self.safety_factor,
                                                                  self.delta_option, self.creator_destructor, self.dem_fem_search,
                                                                  self.search_strategy, self.solver_settings)

        elif self.DEM_parameters["TranslationalIntegrationScheme"].GetString() in {'Hybrid_Bashforth', 'TerminalVelocityScheme'}:
            self.cplusplus_strategy = SDEM.AdamsBashforthStrategy(self.settings, self.max_delta_time, self.n_step_search, self.safety_factor,
                                                                  self.delta_option, self.creator_destructor, self.dem_fem_search,
                                                                  self.search_strategy, self.solver_settings)

        else:
            self.cplusplus_strategy = DEM.ExplicitSolverStrategy(self.settings, self.max_delta_time, self.n_step_search, self.safety_factor,
                                                                 self.delta_option, self.creator_destructor, self.dem_fem_search,
                                                                 self.search_strategy, self.solver_settings)

    def GetTranslationalSchemeInstance(self, class_name):
        try:
            translational_scheme = super(SwimmingStrategy, self).GetTranslationalSchemeInstance(class_name)
        except Exception:
            translational_scheme = SwimmingStrategy.SDEMEvaluateString(class_name)()
        return translational_scheme

    def GetRotationalSchemeInstance(self, class_name):
        try:
            rotational_scheme = super(SwimmingStrategy, self).GetRotationalSchemeInstance(class_name)
        except Exception:
            rotational_scheme = SwimmingStrategy.SDEMEvaluateString(class_name)()
        return rotational_scheme

    def GetHydrodynamicLawParametersIfItExists(self, properties):
        if self.project_parameters.Has('properties'):
            for p in self.project_parameters["properties"]:
                if p['properties_id'].GetInt() == int(properties.Id) and p.Has('hydrodynamic_law_parameters'):
                    return p['hydrodynamic_law_parameters']
        return None

    @staticmethod
    def CreateHydrodynamicLaw(properties, hydrodynamic_law_parameters):

        hydrodynamic_name = hydrodynamic_law_parameters['name'].GetString()
        HydrodynamicInteractionLaw = SwimmingStrategy.SDEMEvaluateString(hydrodynamic_name)(properties, hydrodynamic_law_parameters)

        if hydrodynamic_law_parameters.Has('buoyancy_parameters'):
            buoyancy_parameters = hydrodynamic_law_parameters['buoyancy_parameters']
            buoyancy_name = buoyancy_parameters['name'].GetString()
            if not buoyancy_name == 'default':
                buoyancy_law = SwimmingStrategy.SDEMEvaluateString(buoyancy_name)(buoyancy_parameters)
                HydrodynamicInteractionLaw.SetBuoyancyLaw(buoyancy_law)

        if hydrodynamic_law_parameters.Has('inviscid_force_parameters'):
            inviscid_force_parameters = hydrodynamic_law_parameters['inviscid_force_parameters']
            inviscid_force_name = inviscid_force_parameters['name'].GetString()
            if not inviscid_force_name == 'default':
                inviscid_force_law = SwimmingStrategy.SDEMEvaluateString(inviscid_force_name)(inviscid_force_parameters)
                HydrodynamicInteractionLaw.SetInviscidForceLaw(inviscid_force_law)

        if hydrodynamic_law_parameters.Has('drag_parameters'):
            drag_parameters = hydrodynamic_law_parameters['drag_parameters']
            drag_name = drag_parameters['name'].GetString()
            if not drag_name == 'default':
                drag_law = SwimmingStrategy.SDEMEvaluateString(drag_name)(drag_parameters)
                HydrodynamicInteractionLaw.SetDragLaw(drag_law)

        if hydrodynamic_law_parameters.Has('history_force_parameters'):
            history_force_parameters = hydrodynamic_law_parameters['history_force_parameters']
            history_force_name = history_force_parameters['name'].GetString()
            if not history_force_name == 'default':
                history_force_law = SwimmingStrategy.SDEMEvaluateString(history_force_name)(history_force_parameters)
                HydrodynamicInteractionLaw.SetHistoryForceLaw(history_force_law)

        if hydrodynamic_law_parameters.Has('vorticity_induced_lift_parameters'):
            vorticity_induced_lift_parameters = hydrodynamic_law_parameters['vorticity_induced_lift_parameters']
            vorticity_induced_lift_name = vorticity_induced_lift_parameters['name'].GetString()
            if not vorticity_induced_lift_name == 'default':
                vorticity_induced_lift_law = SwimmingStrategy.SDEMEvaluateString(vorticity_induced_lift_name)(vorticity_induced_lift_parameters)
                HydrodynamicInteractionLaw.SetVorticityInducedLiftLaw(vorticity_induced_lift_law)

        if hydrodynamic_law_parameters.Has('rotation_induced_lift_parameters'):
            rotation_induced_lift_parameters = hydrodynamic_law_parameters['rotation_induced_lift_parameters']
            rotation_induced_lift_name = rotation_induced_lift_parameters['name'].GetString()
            if not rotation_induced_lift_name == 'default':
                rotation_induced_lift_law = SwimmingStrategy.SDEMEvaluateString(rotation_induced_lift_name)(rotation_induced_lift_parameters)
                HydrodynamicInteractionLaw.SetRotationInducedLiftLaw(rotation_induced_lift_law)

        if hydrodynamic_law_parameters.Has('steady_viscous_torque_parameters'):
            steady_viscous_torque_parameters = hydrodynamic_law_parameters['steady_viscous_torque_parameters']
            steady_viscous_torque_name = steady_viscous_torque_parameters['name'].GetString()
            if not steady_viscous_torque_name == 'default':
                steady_viscous_torque_law = SwimmingStrategy.SDEMEvaluateString(steady_viscous_torque_name)(steady_viscous_torque_parameters)
                HydrodynamicInteractionLaw.SetSteadyViscousTorqueLaw(steady_viscous_torque_law)

        HydrodynamicInteractionLaw.SetHydrodynamicInteractionLawInProperties(properties, True)

    def ModifyProperties(self, properties, param = 0):

        super(SwimmingStrategy, self).ModifyProperties(properties, param)

        hydrodynamic_law_parameters = self.GetHydrodynamicLawParametersIfItExists(properties)
        if hydrodynamic_law_parameters:
            SwimmingStrategy.CreateHydrodynamicLaw(properties, hydrodynamic_law_parameters)

        if not param:
            if not properties.Has(Kratos.PARTICLE_SPHERICITY):
                properties[Kratos.PARTICLE_SPHERICITY] = 1.0