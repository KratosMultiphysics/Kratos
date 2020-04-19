from __future__ import print_function, absolute_import, division

# import kratos
import KratosMultiphysics as Kratos

# import RANS
import KratosMultiphysics.RANSApplication as KratosRANS

# import formulation interface
from KratosMultiphysics.RANSApplication.formulations.formulation import Formulation

# import formulations
from ..incompressible_potential_flow import IncompressiblePotentialFlowFormulation
from ..k_epsilon_high_re import KEpsilonHighReFormulation
from .fractional_step_velocity_pressure_formulation import FractionalStepVelocityPressureFormulation

# import utilities
from KratosMultiphysics.RANSApplication import RansVariableUtilities
from KratosMultiphysics.RANSApplication.formulations.utilities import GetConvergenceInfo

class FractionalStepKEpsilonHighReFormulation(Formulation):
    def __init__(self, model_part, settings):
        super(FractionalStepKEpsilonHighReFormulation, self).__init__(model_part, settings)

        default_settings = Kratos.Parameters(r'''
        {
            "formulation_name": "fractional_step_k_epsilon_high_re",
            "incompressible_potential_flow_initialization_settings": {},
            "fractional_step_flow_solver_settings": {},
            "k_epsilon_high_re_solver_settings": {},
            "steady_convergence_settings": {
                "velocity_tolerances": {
                    "relative_tolerance": 1e-3,
                    "absolute_tolerance": 1e-5
                },
                "pressure_tolerances": {
                    "relative_tolerance": 1e-3,
                    "absolute_tolerance": 1e-5
                },
                "turbulent_kinetic_energy_tolerances": {
                    "relative_tolerance": 1e-3,
                    "absolute_tolerance": 1e-5
                },
                "turbulent_energy_dissipation_rate_tolerances": {
                    "relative_tolerance": 1e-3,
                    "absolute_tolerance": 1e-5
                },
                "turbulent_viscosity_tolerances":
                {
                    "relative_tolerance": 1e-3,
                    "absolute_tolerance": 1e-5
                }
            },
            "coupling_settings":
            {
                "velocity_pressure_coupling_iterations": 5,
                "k_epsilon_coupling_iterations": 10
            }
        }''')
        self.settings.ValidateAndAssignDefaults(default_settings)

        if (not self.settings["incompressible_potential_flow_initialization_settings"].IsEquivalentTo(
                Kratos.Parameters("{}"))):
            self.incompressible_potential_flow_formulation = IncompressiblePotentialFlowFormulation(model_part, settings["incompressible_potential_flow_initialization_settings"])
            self.AddFormulation(self.incompressible_potential_flow_formulation)
            self.is_potential_flow_initialized = True

        self.fractional_step_formulation = FractionalStepVelocityPressureFormulation(model_part, settings["fractional_step_flow_solver_settings"])
        self.fractional_step_formulation.SetMaxCouplingIterations(self.settings["coupling_settings"]["velocity_pressure_coupling_iterations"].GetInt())
        self.AddFormulation(self.fractional_step_formulation)

        self.k_epsilon_formulation = KEpsilonHighReFormulation(model_part, settings["k_epsilon_high_re_solver_settings"])
        self.k_epsilon_formulation.SetMaxCouplingIterations(self.settings["coupling_settings"]["k_epsilon_coupling_iterations"].GetInt())
        self.AddFormulation(self.k_epsilon_formulation)

    def IsConverged(self):
        self.is_converged = True

        if (self.is_steady_simulation):
            settings = self.settings["steady_convergence_settings"]
            formulation_converged = self.__CheckTransientConvergence(Kratos.VELOCITY, settings["velocity_tolerances"])
            self.is_converged = self.is_converged and formulation_converged

            formulation_converged = self.__CheckTransientConvergence(Kratos.PRESSURE, settings["pressure_tolerances"])
            self.is_converged = self.is_converged and formulation_converged

            formulation_converged = self.__CheckTransientConvergence(KratosRANS.TURBULENT_KINETIC_ENERGY, settings["turbulent_kinetic_energy_tolerances"])
            self.is_converged = self.is_converged and formulation_converged

            formulation_converged = self.__CheckTransientConvergence(KratosRANS.TURBULENT_ENERGY_DISSIPATION_RATE, settings["turbulent_energy_dissipation_rate_tolerances"])
            self.is_converged = self.is_converged and formulation_converged

            formulation_converged = self.__CheckTransientConvergence(Kratos.TURBULENT_VISCOSITY, settings["turbulent_viscosity_tolerances"])
            self.is_converged = self.is_converged and formulation_converged
        else:
            self.is_converged = self.k_epsilon_formulation.IsConverged()

        return self.is_converged

    def SetConstants(self, settings):
        self.k_epsilon_formulation.SetConstants(settings)

    def SetTimeSchemeSettings(self, settings):
        scheme_type = settings["scheme_type"].GetString()
        if (scheme_type == "steady"):
            self.is_steady_simulation = True
            default_settings = Kratos.Parameters('''{
                "scheme_type": "steady",
                "pressure_coefficient": 0.5
            }''')
            settings.ValidateAndAssignDefaults(default_settings)
            self.GetBaseModelPart().ProcessInfo.SetValue(Kratos.PRESSURE_COEFFICIENT, settings["pressure_coefficient"].GetDouble())
        elif (scheme_type == "transient"):
            self.is_steady_simulation = False
            default_settings = Kratos.Parameters('''{
                "scheme_type": "transient",
                "alpha_bossak": -0.3
            }''')
            settings.ValidateAndAssignDefaults(default_settings)
            self.GetBaseModelPart().ProcessInfo.SetValue(Kratos.BOSSAK_ALPHA, settings["alpha_bossak"].GetDouble())
        else:
            raise Exception("Only \"steady\" and \"transient\" time schemes are supported by " + self.GetName() + " [ scheme_type = " + scheme_type + " ]")

        super(FractionalStepKEpsilonHighReFormulation, self).SetTimeSchemeSettings(settings)

    def __CheckTransientConvergence(self, variable, settings):
        relative_error, absolute_error = RansVariableUtilities.CalculateTransientVariableConvergence(self.GetBaseModelPart(), variable)
        relative_tolerance = settings["relative_tolerance"].GetDouble()
        absolute_tolerance = settings["absolute_tolerance"].GetDouble()

        info = GetConvergenceInfo(variable, relative_error, relative_tolerance, absolute_error, absolute_tolerance)
        Kratos.Logger.PrintInfo(self.GetName(), info)

        return (relative_error <= relative_tolerance or absolute_error <= absolute_tolerance)