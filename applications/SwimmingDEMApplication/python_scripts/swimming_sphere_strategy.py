import KratosMultiphysics as Kratos
import KratosMultiphysics.DEMApplication as DEM
import KratosMultiphysics.SwimmingDEMApplication as SDEM

from KratosMultiphysics.DEMApplication.sphere_strategy import ExplicitStrategy
BaseStrategy = ExplicitStrategy

def GetTerminalVelocitySchemeParameters(sdem_parameters):
    """Parameters of the TerminalVelocityScheme read from the SwimmingDEM ProjectParameters.

    The scheme moves inertia-free particles with the fluid velocity plus the Stokes
    settling velocity 2 a^2 (rho_p - rho_f) g / (9 mu). Both mu and g must be given
    explicitly in ``custom_dem.terminal_velocity_scheme_parameters``::

        "custom_dem" : {
            "translational_integration_scheme" : "TerminalVelocityScheme",
            "terminal_velocity_scheme_parameters" : {
                "dynamic_viscosity" : 1.0,
                "gravity"           : [0.0, 0.0, -1.0]
            }
        }

    in the units of the case (dimensional or dimensionless); the scheme itself
    carries no hidden constants. The gravity is checked against ``gravity_parameters``
    (the gravity seen by the fluid) so that both cannot silently disagree.
    """
    if not sdem_parameters.Has("custom_dem") or not sdem_parameters["custom_dem"].Has("terminal_velocity_scheme_parameters"):
        raise Exception("TerminalVelocityScheme requires 'custom_dem.terminal_velocity_scheme_parameters' "
                        "with 'dynamic_viscosity' and 'gravity' in the SwimmingDEM ProjectParameters "
                        "(no hidden default values are used).")
    scheme_parameters = sdem_parameters["custom_dem"]["terminal_velocity_scheme_parameters"]
    for key in ("dynamic_viscosity", "gravity"):
        if not scheme_parameters.Has(key):
            raise Exception("'custom_dem.terminal_velocity_scheme_parameters' must contain '{}'".format(key))
    if sdem_parameters.Has("gravity_parameters"):
        g_fluid = sdem_parameters["gravity_parameters"]["direction"].GetVector()
        g_fluid *= sdem_parameters["gravity_parameters"]["modulus"].GetDouble()
        g_scheme = scheme_parameters["gravity"].GetVector()
        scale = max(1.0e-300, max(abs(g_fluid[k]) for k in range(3)), max(abs(g_scheme[k]) for k in range(3)))
        if max(abs(g_fluid[k] - g_scheme[k]) for k in range(3)) > 1.0e-10 * scale:
            Kratos.Logger.PrintWarning("TerminalVelocityScheme",
                "the gravity of 'terminal_velocity_scheme_parameters' {} differs from 'gravity_parameters' {}".format(
                    [g_scheme[k] for k in range(3)], [g_fluid[k] for k in range(3)]))
    return scheme_parameters


def CreateTerminalVelocityScheme(sdem_parameters):
    return SDEM.TerminalVelocityScheme(GetTerminalVelocitySchemeParameters(sdem_parameters))


class SwimmingStrategy(BaseStrategy):

    @staticmethod
    def SDEMEvaluateString(name):
        return getattr(SDEM, name)

    def __init__(self, all_model_parts, creator_destructor, dem_fem_search, parameters, procedures):
        self.project_parameters = parameters
        super().__init__(all_model_parts, creator_destructor, dem_fem_search, parameters['dem_parameters'], procedures)

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

        if self.DEM_parameters["TranslationalIntegrationScheme"].GetString() == 'Velocity_Verlet':
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
        if class_name == 'TerminalVelocityScheme':
            return CreateTerminalVelocityScheme(self.project_parameters)
        try:
            translational_scheme = super().GetTranslationalSchemeInstance(class_name)
        except Exception:
            translational_scheme = SwimmingStrategy.SDEMEvaluateString(class_name)()
        return translational_scheme

    def GetRotationalSchemeInstance(self, class_name):
        if class_name == 'TerminalVelocityScheme':
            return CreateTerminalVelocityScheme(self.project_parameters)
        try:
            rotational_scheme = super().GetRotationalSchemeInstance(class_name)
        except Exception:
            rotational_scheme = SwimmingStrategy.SDEMEvaluateString(class_name)()
        return rotational_scheme

    def GetHydrodynamicLawParametersIfItExists(self, properties):
        if self.project_parameters.Has('properties'):
            for p in self.project_parameters["properties"]:
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

        if hydrodynamic_law_parameters.Has('virtual_mass_force_parameters'):
            virtual_mass_force_parameters = hydrodynamic_law_parameters['virtual_mass_force_parameters']
            virtual_mass_force_name = virtual_mass_force_parameters['name'].GetString()
            if not virtual_mass_force_name == 'default':
                virtual_mass_force_law = SwimmingStrategy.SDEMEvaluateString(virtual_mass_force_name)(virtual_mass_force_parameters)
                HydrodynamicInteractionLaw.SetVirtualMassForceLaw(virtual_mass_force_law)

        if hydrodynamic_law_parameters.Has('undisturbed_force_parameters'):
            undisturbed_force_parameters = hydrodynamic_law_parameters['undisturbed_force_parameters']
            undisturbed_force_name = undisturbed_force_parameters['name'].GetString()
            if not undisturbed_force_name == 'default':
                undisturbed_force_law = SwimmingStrategy.SDEMEvaluateString(undisturbed_force_name)(undisturbed_force_parameters)
                HydrodynamicInteractionLaw.SetUndisturbedForceLaw(undisturbed_force_law)

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

        super().ModifyProperties(properties, param)

        hydrodynamic_law_parameters = self.GetHydrodynamicLawParametersIfItExists(properties)
        if hydrodynamic_law_parameters:
            SwimmingStrategy.CreateHydrodynamicLaw(properties, hydrodynamic_law_parameters)

        if not param:
            if not properties.Has(Kratos.PARTICLE_SPHERICITY):
                properties[Kratos.PARTICLE_SPHERICITY] = 1.0