from __future__ import print_function, absolute_import, division

# import kratos
import KratosMultiphysics as Kratos

# import required applications
import KratosMultiphysics.RANSApplication as KratosRANS

# import formulation interface
from KratosMultiphysics.RANSApplication.formulations.rans_formulation import RansFormulation

# import utilities
from KratosMultiphysics.RANSApplication.formulations.utilities import CreateBlockBuilderAndSolver
from KratosMultiphysics.RANSApplication.formulations.utilities import CreateRansFormulationModelPart
from KratosMultiphysics.RANSApplication.formulations.utilities import CreateAlgebraicFluxCorrectedSteadyScalarScheme
from KratosMultiphysics.RANSApplication.formulations.utilities import CalculateNormalsOnConditions
from KratosMultiphysics.RANSApplication.formulations.utilities import InitializePeriodicConditions
from KratosMultiphysics.RANSApplication.formulations.utilities import InitializeYPlusVariablesInConditions
from KratosMultiphysics.RANSApplication.formulations.utilities import GetKratosObjectType
from KratosMultiphysics.RANSApplication.formulations.utilities import GetBoundaryFlags

class KOmegaOmegaRansFormulation(RansFormulation):
    def __init__(self, model_part, settings):
        super().__init__(model_part, settings)

        defaults = Kratos.Parameters(r"""{
            "relative_tolerance"    : 1e-3,
            "absolute_tolerance"    : 1e-5,
            "max_iterations"        : 200,
            "relaxation_factor"     : 0.5,
            "echo_level"            : 2,
            "linear_solver_settings": {
                "solver_type"  : "amgcl"
            },
            "boundary_flags": ["INLET", "STRUCTURE"]
        }""")

        settings = self.GetParameters()
        settings.ValidateAndAssignDefaults(defaults)
        self.echo_level = settings["echo_level"].GetInt()

    def PrepareModelPart(self):
        self.omega_model_part = CreateRansFormulationModelPart(self, self.element_name, self.condition_name)

        Kratos.Logger.PrintInfo(self.GetName(),
                                "Created formulation model part.")

    def Initialize(self):
        InitializeYPlusVariablesInConditions(self.omega_model_part)
        CalculateNormalsOnConditions(self.omega_model_part)

        if (self.IsPeriodic()):
            InitializePeriodicConditions(self.GetBaseModelPart(), self.omega_model_part, [KratosRANS.TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE])

        solver_settings = self.GetParameters()

        linear_solver = GetKratosObjectType("LinearSolverFactory")(
            solver_settings["linear_solver_settings"])

        builder_and_solver = CreateBlockBuilderAndSolver(
            linear_solver, self.IsPeriodic(), self.GetCommunicator())

        convergence_criteria = GetKratosObjectType("MixedGenericCriteria")([
            (KratosRANS.TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE,
            solver_settings["relative_tolerance"].GetDouble(),
            solver_settings["absolute_tolerance"].GetDouble())
        ])

        if (self.is_steady_simulation):
            scheme = self.scheme_type(solver_settings["relaxation_factor"].GetDouble())
        else:
            scheme = GetKratosObjectType("BossakRelaxationScalarScheme")(
                self.omega_model_part.ProcessInfo[Kratos.BOSSAK_ALPHA],
                solver_settings["relaxation_factor"].GetDouble(),
                KratosRANS.TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE)

        self.solver = GetKratosObjectType("ResidualBasedNewtonRaphsonStrategy")(
            self.omega_model_part, scheme,
            convergence_criteria, builder_and_solver, solver_settings["max_iterations"].GetInt(), False,
            False, False)

        builder_and_solver.SetEchoLevel(
            solver_settings["echo_level"].GetInt() - 3)
        self.solver.SetEchoLevel(
            solver_settings["echo_level"].GetInt() - 2)
        convergence_criteria.SetEchoLevel(
            solver_settings["echo_level"].GetInt() - 1)

        super().Initialize()
        Kratos.Logger.PrintInfo(self.GetName(), "Initialized formulation")

    def SolveCouplingStep(self):
        if (self.IsBufferInitialized()):
            self.solver.Predict()
            self.solver.SolveSolutionStep()
            Kratos.Logger.PrintInfo(self.GetName(), "Solved  formulation.")
            return True

        return False

    def GetStrategy(self):
        return self.solver

    def SetTimeSchemeSettings(self, settings):
        if (settings.Has("scheme_type")):
            scheme_type = settings["scheme_type"].GetString()
            if (scheme_type == "steady"):
                self.is_steady_simulation = True
            elif (scheme_type == "bdf2" or scheme_type == "bossak"):
                self.is_steady_simulation = False
            else:
                raise Exception(
                    "Only \"steady\", \"bdf2\" and \"bossak\" scheme types supported. [ scheme_type = \""
                    + scheme_type + "\" ]")
        else:
            raise Exception(
                "\"scheme_type\" is missing in time scheme settings")

    def GetMaxCouplingIterations(self):
        return "N/A"

    def GetModelPart(self):
        return self.omega_model_part

    def SetStabilizationMethod(self, stabilization_method):
        if (stabilization_method == "algebraic_flux_corrected"):
            self.element_name = "RansKOmegaOmegaAFC"
            self.scheme_type = lambda x: CreateAlgebraicFluxCorrectedSteadyScalarScheme(
                                                    x,
                                                    GetBoundaryFlags(self.GetParameters()["boundary_flags"]),
                                                    self.IsPeriodic())
        elif (stabilization_method == "residual_based_flux_corrected"):
            self.element_name = "RansKOmegaOmegaRFC"
            self.scheme_type = GetKratosObjectType("SteadyScalarScheme")
        elif (stabilization_method == "non_linear_cross_wind_dissipation"):
            self.element_name = "RansKOmegaOmegaCWD"
            self.scheme_type = GetKratosObjectType("SteadyScalarScheme")
        else:
            raise Exception("Unsupported stabilization method")

    def SetWallFunctionSettings(self, settings):
        wall_function_region_type = "logarithmic_region_only"
        if (settings.Has("wall_function_region_type")):
            wall_function_region_type = settings[
                "wall_function_region_type"].GetString()

        wall_friction_velocity_calculation_method = "velocity_based"
        if (settings.Has("wall_friction_velocity_calculation_method")):
            wall_friction_velocity_calculation_method = settings[
                "wall_friction_velocity_calculation_method"].GetString()

        if (wall_function_region_type == "logarithmic_region_only"):
            if (wall_friction_velocity_calculation_method == "velocity_based"):
                self.condition_name = "RansKOmegaOmegaUBasedWall"
            elif (wall_friction_velocity_calculation_method ==
                  "turbulent_kinetic_energy_based"):
                self.condition_name = "RansKOmegaOmegaKBasedWall"
            else:
                msg = "Unsupported wall friction velocity calculation method. [ wall_friction_velocity_calculation_method = \"" + wall_friction_velocity_calculation_method + "\" ].\n"
                msg += "Supported methods are:\n"
                msg += "\tvelocity_based\n"
                msg += "\tturbulent_kinetic_energy_based\n"
                raise Exception(msg)
        else:
            msg = "Unsupported wall function region type provided. [ wall_function_region_type = \"" + wall_function_region_type + "\" ]."
            msg += "Supported wall function region types are:\n"
            msg += "\tlogarithmic_region_only\n"
            raise Exception(msg)