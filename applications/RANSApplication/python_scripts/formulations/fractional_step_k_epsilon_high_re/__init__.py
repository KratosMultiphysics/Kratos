from __future__ import print_function, absolute_import, division

# import kratos
import KratosMultiphysics as Kratos

# import RANS
import KratosMultiphysics.FluidDynamicsApplication as KratosCFD
import KratosMultiphysics.RANSApplication as KratosRANS

# import formulation interface
from KratosMultiphysics.RANSApplication.formulations.formulation import Formulation

# import formulations
from ..incompressible_potential_flow import IncompressiblePotentialFlowFormulation
from ..fractional_step_velocity_pressure_formulation import FractionalStepVelocityPressureFormulation
from .k_epsilon_high_re_k_formulation import KEpsilonHighReKFormulation
from .k_epsilon_high_re_epsilon_formulation import KEpsilonHighReEpsilonFormulation

# import utilities
from KratosMultiphysics.RANSApplication import RansCalculationUtilities
from KratosMultiphysics.RANSApplication import RansVariableUtilities
from KratosMultiphysics.RANSApplication.formulations.utilities import GetConvergenceInfo

class FractionalStepKEpsilonhighReFormulation(Formulation):
    def __init__(self, model_part, settings):
        super(FractionalStepKEpsilonhighReFormulation, self).__init__(model_part, settings)

        default_settings = Kratos.Parameters(r'''
        {
            "formulation_name": "fractional_step_k_epsilon_high_re",
            "incompressible_potential_flow_initialization_settings": {},
            "fractional_step_flow_solver_settings": {},
            "turbulent_kinetic_energy_solver_settings": {},
            "turbulent_energy_dissipation_rate_solver_settings": {},
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
        self.AddFormulation(self.fractional_step_formulation)

        self.tke_formulation = KEpsilonHighReKFormulation(model_part, settings["turbulent_kinetic_energy_solver_settings"])
        self.AddFormulation(self.tke_formulation)

        self.epsilon_formulation = KEpsilonHighReEpsilonFormulation(model_part, settings["turbulent_energy_dissipation_rate_solver_settings"])
        self.AddFormulation(self.epsilon_formulation)

    def AddVariables(self):
        self.GetBaseModelPart().AddNodalSolutionStepVariable(Kratos.DENSITY)
        self.GetBaseModelPart().AddNodalSolutionStepVariable(Kratos.PRESSURE)
        self.GetBaseModelPart().AddNodalSolutionStepVariable(Kratos.VELOCITY)
        self.GetBaseModelPart().AddNodalSolutionStepVariable(Kratos.ACCELERATION)
        self.GetBaseModelPart().AddNodalSolutionStepVariable(Kratos.MESH_VELOCITY)
        self.GetBaseModelPart().AddNodalSolutionStepVariable(Kratos.BODY_FORCE)
        self.GetBaseModelPart().AddNodalSolutionStepVariable(Kratos.NODAL_H)
        self.GetBaseModelPart().AddNodalSolutionStepVariable(Kratos.REACTION)
        self.GetBaseModelPart().AddNodalSolutionStepVariable(Kratos.REACTION_WATER_PRESSURE)
        self.GetBaseModelPart().AddNodalSolutionStepVariable(Kratos.NORMAL)
        self.GetBaseModelPart().AddNodalSolutionStepVariable(Kratos.Y_WALL)
        self.GetBaseModelPart().AddNodalSolutionStepVariable(Kratos.EXTERNAL_PRESSURE)
        self.GetBaseModelPart().AddNodalSolutionStepVariable(Kratos.VISCOSITY)
        self.GetBaseModelPart().AddNodalSolutionStepVariable(Kratos.KINEMATIC_VISCOSITY)
        self.GetBaseModelPart().AddNodalSolutionStepVariable(Kratos.FRACT_VEL)
        self.GetBaseModelPart().AddNodalSolutionStepVariable(Kratos.PRESSURE_OLD_IT)
        # The following are used for the calculation of projections
        self.GetBaseModelPart().AddNodalSolutionStepVariable(Kratos.NODAL_AREA)
        self.GetBaseModelPart().AddNodalSolutionStepVariable(Kratos.PRESS_PROJ)
        self.GetBaseModelPart().AddNodalSolutionStepVariable(Kratos.CONV_PROJ)
        self.GetBaseModelPart().AddNodalSolutionStepVariable(Kratos.DIVPROJ)
        self.GetBaseModelPart().AddNodalSolutionStepVariable(KratosCFD.Q_VALUE)

        self.GetBaseModelPart().AddNodalSolutionStepVariable(Kratos.TURBULENT_VISCOSITY)
        self.GetBaseModelPart().AddNodalSolutionStepVariable(KratosRANS.RANS_Y_PLUS)
        self.GetBaseModelPart().AddNodalSolutionStepVariable(KratosRANS.TURBULENT_KINETIC_ENERGY)
        self.GetBaseModelPart().AddNodalSolutionStepVariable(KratosRANS.TURBULENT_KINETIC_ENERGY_RATE)
        self.GetBaseModelPart().AddNodalSolutionStepVariable(KratosRANS.TURBULENT_ENERGY_DISSIPATION_RATE)
        self.GetBaseModelPart().AddNodalSolutionStepVariable(KratosRANS.TURBULENT_ENERGY_DISSIPATION_RATE_2)
        self.GetBaseModelPart().AddNodalSolutionStepVariable(KratosRANS.RANS_AUXILIARY_VARIABLE_1)
        self.GetBaseModelPart().AddNodalSolutionStepVariable(KratosRANS.RANS_AUXILIARY_VARIABLE_2)

        if (hasattr(self, "is_potential_flow_initialized")):
            self.GetBaseModelPart().AddNodalSolutionStepVariable(KratosRANS.VELOCITY_POTENTIAL)
            self.GetBaseModelPart().AddNodalSolutionStepVariable(KratosRANS.PRESSURE_POTENTIAL)

        Kratos.Logger.PrintInfo(self.GetName(), "Added solution step variables.")

    def AddDofs(self):
        Kratos.VariableUtils().AddDof(Kratos.VELOCITY_X, Kratos.REACTION_X,self.GetBaseModelPart())
        Kratos.VariableUtils().AddDof(Kratos.VELOCITY_Y, Kratos.REACTION_Y,self.GetBaseModelPart())
        Kratos.VariableUtils().AddDof(Kratos.VELOCITY_Z, Kratos.REACTION_Z,self.GetBaseModelPart())
        Kratos.VariableUtils().AddDof(Kratos.PRESSURE, Kratos.REACTION_WATER_PRESSURE,self.GetBaseModelPart())
        Kratos.VariableUtils().AddDof(KratosRANS.TURBULENT_KINETIC_ENERGY, self.GetBaseModelPart())
        Kratos.VariableUtils().AddDof(KratosRANS.TURBULENT_ENERGY_DISSIPATION_RATE, self.GetBaseModelPart())
        if (hasattr(self, "is_potential_flow_initialized")):
            Kratos.VariableUtils().AddDof(KratosRANS.VELOCITY_POTENTIAL, self.GetBaseModelPart())
            Kratos.VariableUtils().AddDof(KratosRANS.PRESSURE_POTENTIAL, self.GetBaseModelPart())

        Kratos.Logger.PrintInfo(self.GetName(), "Added solution step dofs.")

    def GetMinimumBufferSize(self):
        return 3

    def IsConverged(self):
        self.is_converged = True

        strategy_converged = self.fractional_step_formulation.GetStrategy().IsConverged()
        if (self.is_converged):
            self.is_converged = strategy_converged

        strategy_converged = self.tke_formulation.GetStrategy().IsConverged()
        if (self.is_converged):
            self.is_converged = strategy_converged

        strategy_converged = self.epsilon_formulation.GetStrategy().IsConverged()
        if (self.is_converged):
            self.is_converged = strategy_converged

        if (self.is_steady_simulation):
            settings = self.settings["steady_convergence_settings"]
            formulation_converged = self.__CheckTransientConvergence(Kratos.VELOCITY, settings["velocity_tolerances"])
            if (self.is_converged):
                self.is_converged = formulation_converged

            formulation_converged = self.__CheckTransientConvergence(Kratos.PRESSURE, settings["pressure_tolerances"])
            if (self.is_converged):
                self.is_converged = formulation_converged

            formulation_converged = self.__CheckTransientConvergence(KratosRANS.TURBULENT_KINETIC_ENERGY, settings["turbulent_kinetic_energy_tolerances"])
            if (self.is_converged):
                self.is_converged = formulation_converged

            formulation_converged = self.__CheckTransientConvergence(KratosRANS.TURBULENT_ENERGY_DISSIPATION_RATE, settings["turbulent_energy_dissipation_rate_tolerances"])
            if (self.is_converged):
                self.is_converged = formulation_converged

        return self.is_converged


    def SolveCouplingStep(self):
        if (not hasattr(self, "is_initialized")):
            self.is_initialized = True
            if (hasattr(self, "incompressible_potential_flow_formulation")):
                self.incompressible_potential_flow_formulation.SolveCouplingStep()
                Kratos.Logger.PrintInfo(self.GetName(), "Initialized with incompressible potential flow solver")

        settings = self.settings["coupling_settings"]
        max_iterations = settings["velocity_pressure_coupling_iterations"].GetInt()
        for iteration in range(max_iterations):
            self.fractional_step_formulation.ExecuteBeforeCouplingSolveStep()
            self.fractional_step_formulation.SolveCouplingStep()
            self.fractional_step_formulation.ExecuteAfterCouplingSolveStep()
            Kratos.Logger.PrintInfo(self.GetName(), "Solved itr. " + str(iteration + 1) + "/" + str(max_iterations) + " velocity-pressure coupling.")

        max_iterations = settings["k_epsilon_coupling_iterations"].GetInt()
        for iteration in range(max_iterations):
            self.tke_formulation.ExecuteBeforeCouplingSolveStep()
            self.tke_formulation.SolveCouplingStep()
            self.tke_formulation.ExecuteAfterCouplingSolveStep()
            Kratos.Logger.PrintInfo(self.GetName(), "Solved k equation.")
            self.epsilon_formulation.ExecuteBeforeCouplingSolveStep()
            self.epsilon_formulation.SolveCouplingStep()
            self.epsilon_formulation.ExecuteAfterCouplingSolveStep()
            Kratos.Logger.PrintInfo(self.GetName(), "Solved epsilon equation.")
            Kratos.Logger.PrintInfo(self.GetName(), "Solved itr. " + str(iteration + 1) + "/" + str(max_iterations) + " k-epsilon coupling.")

    def SetConstants(self, settings):
        defaults = Kratos.Parameters('''{
            "wall_smoothness_beta"    : 5.2,
            "von_karman"              : 0.41,
            "c_mu"                    : 0.09,
            "c1"                      : 1.44,
            "c2"                      : 1.92,
            "sigma_k"                 : 1.0,
            "sigma_epsilon"           : 1.3,
            "stabilizing_upwind_operator_coefficient": 1.2,
            "stabilizing_positivity_preserving_coefficient": 1.2
        }''')

        settings.ValidateAndAssignDefaults(defaults)

        process_info = self.GetBaseModelPart().ProcessInfo
        process_info.SetValue(KratosRANS.WALL_SMOOTHNESS_BETA, settings["wall_smoothness_beta"].GetDouble())
        process_info.SetValue(KratosRANS.WALL_VON_KARMAN, settings["von_karman"].GetDouble())
        process_info.SetValue(KratosRANS.TURBULENCE_RANS_C_MU, settings["c_mu"].GetDouble())
        process_info.SetValue(KratosRANS.TURBULENCE_RANS_C1, settings["c1"].GetDouble())
        process_info.SetValue(KratosRANS.TURBULENCE_RANS_C2, settings["c2"].GetDouble())
        process_info.SetValue(KratosRANS.TURBULENT_KINETIC_ENERGY_SIGMA, settings["sigma_k"].GetDouble())
        process_info.SetValue(KratosRANS.TURBULENT_ENERGY_DISSIPATION_RATE_SIGMA, settings["sigma_epsilon"].GetDouble())
        process_info.SetValue(KratosRANS.RANS_STABILIZATION_DISCRETE_UPWIND_OPERATOR_COEFFICIENT, settings["stabilizing_upwind_operator_coefficient"].GetDouble())
        process_info.SetValue(KratosRANS.RANS_STABILIZATION_DIAGONAL_POSITIVITY_PRESERVING_COEFFICIENT, settings["stabilizing_positivity_preserving_coefficient"].GetDouble())
        process_info.SetValue(KratosRANS.RANS_Y_PLUS_LIMIT, RansCalculationUtilities.CalculateLogarithmicYPlusLimit(
                                                                                    process_info[KratosRANS.WALL_VON_KARMAN],
                                                                                    process_info[KratosRANS.WALL_SMOOTHNESS_BETA]
                                                                                    ))

    def SetTimeSchemeSettings(self, settings):
        scheme_type = settings["scheme_type"].GetString()
        if (scheme_type == "steady"):
            self.is_steady_simulation = True
            self.fractional_step_formulation.is_steady_simulation = True
            self.tke_formulation.is_steady_simulation = True
            self.epsilon_formulation.is_steady_simulation = True
            default_settings = Kratos.Parameters('''{
                "scheme_type": "steady",
                "pressure_coefficient": 0.5
            }''')
            settings.ValidateAndAssignDefaults(default_settings)
            self.GetBaseModelPart().ProcessInfo.SetValue(Kratos.PRESSURE_COEFFICIENT, settings["pressure_coefficient"].GetDouble())
        elif (scheme_type == "transient"):
            self.is_steady_simulation = False
            self.fractional_step_formulation.is_steady_simulation = False
            self.tke_formulation.is_steady_simulation = False
            self.epsilon_formulation.is_steady_simulation = False
            default_settings = Kratos.Parameters('''{
                "scheme_type": "transient",
                "alpha_bossak": -0.3
            }''')
            settings.ValidateAndAssignDefaults(default_settings)
            self.GetBaseModelPart().ProcessInfo.SetValue(Kratos.BOSSAK_ALPHA, settings["alpha_bossak"].GetDouble())
        else:
            raise Exception("Only \"steady\" and \"transient\" time schemes are supported by " + self.GetName() + " [ scheme_type = " + scheme_type + " ]")

        self.time_scheme_settings = settings

    def __CheckTransientConvergence(self, variable, settings):
        relative_error, absolute_error = RansVariableUtilities.CalculateTransientVariableConvergence(self.GetBaseModelPart(), variable)
        relative_tolerance = settings["relative_tolerance"].GetDouble()
        absolute_tolerance = settings["absolute_tolerance"].GetDouble()

        info = GetConvergenceInfo(variable, relative_error, relative_tolerance, absolute_error, absolute_tolerance)
        Kratos.Logger.PrintInfo(self.GetName(), info)

        return (relative_error <= relative_tolerance or absolute_error <= absolute_tolerance)