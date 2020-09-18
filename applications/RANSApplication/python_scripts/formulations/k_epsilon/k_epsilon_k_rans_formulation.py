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
from KratosMultiphysics.RANSApplication.formulations.utilities import InitializePeriodicConditions
from KratosMultiphysics.RANSApplication.formulations.utilities import GetKratosObjectType
from KratosMultiphysics.RANSApplication.formulations.utilities import GetBoundaryFlags

class KEpsilonKRansFormulation(RansFormulation):
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
        self.k_model_part = CreateRansFormulationModelPart(self, self.element_name)

        Kratos.Logger.PrintInfo(self.GetName(),
                                "Created formulation model part.")

    def Initialize(self):
        solver_settings = self.GetParameters()

        if (self.IsPeriodic()):
            InitializePeriodicConditions(self.GetBaseModelPart(), self.k_model_part, [KratosRANS.TURBULENT_KINETIC_ENERGY])

        solver_settings = self.GetParameters()

        linear_solver = GetKratosObjectType("LinearSolverFactory")(
            solver_settings["linear_solver_settings"])

        builder_and_solver = CreateBlockBuilderAndSolver(
            linear_solver, self.IsPeriodic(), self.GetCommunicator())

        convergence_criteria = GetKratosObjectType("MixedGenericCriteria")([
            (KratosRANS.TURBULENT_KINETIC_ENERGY,
            solver_settings["relative_tolerance"].GetDouble(),
            solver_settings["absolute_tolerance"].GetDouble())
        ])

        if (self.is_steady_simulation):
            scheme = self.scheme_type(solver_settings["relaxation_factor"].GetDouble())
        else:
            scheme = GetKratosObjectType("BossakRelaxationScalarScheme")(
                self.k_model_part.ProcessInfo[Kratos.BOSSAK_ALPHA],
                solver_settings["relaxation_factor"].GetDouble(),
                KratosRANS.TURBULENT_KINETIC_ENERGY,
                KratosRANS.TURBULENT_KINETIC_ENERGY_RATE,
                KratosRANS.RANS_AUXILIARY_VARIABLE_1)

        self.solver = GetKratosObjectType("ResidualBasedNewtonRaphsonStrategy")(
            self.k_model_part, scheme,
            convergence_criteria, builder_and_solver, solver_settings["max_iterations"].GetInt(), False,
            False, False)

        builder_and_solver.SetEchoLevel(
            solver_settings["echo_level"].GetInt() - 3)
        self.solver.SetEchoLevel(
            solver_settings["echo_level"].GetInt() - 2)
        convergence_criteria.SetEchoLevel(
            solver_settings["echo_level"].GetInt() - 1)

        self.solver.Initialize()
        Kratos.Logger.PrintInfo(self.GetName(), "Initialized formulation")

    def InitializeSolutionStep(self):
        self.solver.InitializeSolutionStep()

    def SolveCouplingStep(self):
        if (self.IsBufferInitialized()):
            self.solver.Predict()
            self.solver.SolveSolutionStep()
            Kratos.Logger.PrintInfo(self.GetName(), "Solved  formulation.")
            return True

        return False

    def FinalizeSolutionStep(self):
        self.solver.FinalizeSolutionStep()

    def Check(self):
        self.solver.Check()

    def Clear(self):
        self.solver.Clear()

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
            raise Exception("\"scheme_type\" is missing in time scheme settings")

    def GetMaxCouplingIterations(self):
        return "N/A"

    def GetModelPart(self):
        return self.k_model_part

    def SetStabilizationMethod(self, stabilization_method):
        if (stabilization_method == "algebraic_flux_corrected"):
            self.element_name = "RansKEpsilonKAFC"
            self.scheme_type = lambda x: CreateAlgebraicFluxCorrectedSteadyScalarScheme(
                                                    x,
                                                    GetBoundaryFlags(self.GetParameters()["boundary_flags"]),
                                                    self.IsPeriodic())
        elif (stabilization_method == "residual_based_flux_corrected"):
            self.element_name = "RansKEpsilonKRFC"
            self.scheme_type = GetKratosObjectType("SteadyScalarScheme")
        elif (stabilization_method == "non_linear_cross_wind_dissipation"):
            self.element_name = "RansKEpsilonKCWD"
            self.scheme_type = GetKratosObjectType("SteadyScalarScheme")
        else:
            raise Exception("Unsupported stabilization method")