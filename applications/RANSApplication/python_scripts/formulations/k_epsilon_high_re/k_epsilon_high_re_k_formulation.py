from __future__ import print_function, absolute_import, division

# import kratos
import KratosMultiphysics as Kratos

# import required applications
import KratosMultiphysics.RANSApplication as KratosRANS

# import formulation interface
from KratosMultiphysics.RANSApplication.formulations.formulation import Formulation

# import utilities
from KratosMultiphysics.RANSApplication.formulations.utilities import CreateLinearSolver
from KratosMultiphysics.RANSApplication.formulations.utilities import CreateFormulationModelPart
from KratosMultiphysics.RANSApplication.formulations.utilities import CalculateNormalsOnConditions
from KratosMultiphysics.RANSApplication.formulations.utilities import CreateResidualBasedBlockBuilderAndSolver
from KratosMultiphysics.RANSApplication.formulations.utilities import CreateResidualCriteria
from KratosMultiphysics.RANSApplication.formulations.utilities import CreateResidualBasedNewtonRaphsonStrategy
from KratosMultiphysics.RANSApplication.formulations.utilities import CreateSteadyAlgeraicFluxCorrectedTransportScheme
from KratosMultiphysics.RANSApplication.formulations.utilities import CreateBossakScalarScheme
from KratosMultiphysics.RANSApplication.formulations.utilities import CreateSteadyScalarScheme
from KratosMultiphysics.RANSApplication.formulations.utilities import IsBufferInitialized

class KEpsilonHighReKFormulation(Formulation):
    def __init__(self, model_part, settings):
        super(KEpsilonHighReKFormulation,
              self).__init__(model_part, settings)

        defaults = Kratos.Parameters(r"""{
            "relative_tolerance"    : 1e-3,
            "absolute_tolerance"    : 1e-5,
            "max_iterations"        : 200,
            "relaxation_factor"     : 0.5,
            "echo_level"            : 2,
            "linear_solver_settings": {
                "solver_type"  : "amgcl"
            }
        }""")

        self.settings.ValidateAndAssignDefaults(defaults)
        self.echo_level = self.settings["echo_level"].GetInt()

    def PrepareModelPart(self):
        if self.GetBaseModelPart().ProcessInfo[Kratos.DOMAIN_SIZE] == 2:
            condition_name = "LineCondition"
        else:
            condition_name = "SurfaceCondition"

        self.k_model_part = CreateFormulationModelPart(self, self.element_name, condition_name)

        Kratos.Logger.PrintInfo(self.GetName(),
                                "Created formulation model part.")

    def Initialize(self):
        solver_settings = self.settings
        linear_solver = CreateLinearSolver(
            solver_settings["linear_solver_settings"])
        builder_and_solver = CreateResidualBasedBlockBuilderAndSolver(
            linear_solver, self.IsPeriodic(), self.GetCommunicator())
        convergence_criteria = CreateResidualCriteria(
                                self.settings["relative_tolerance"].GetDouble(),
                                self.settings["absolute_tolerance"].GetDouble())

        if (self.is_steady_simulation):
            scheme = self.scheme_type(self.settings["relaxation_factor"].GetDouble())
        else:
            scheme = CreateBossakScalarScheme(
                self.k_model_part.ProcessInfo[Kratos.BOSSAK_ALPHA],
                self.settings["relaxation_factor"].GetDouble(),
                KratosRANS.TURBULENT_KINETIC_ENERGY,
                KratosRANS.TURBULENT_KINETIC_ENERGY_RATE,
                KratosRANS.RANS_AUXILIARY_VARIABLE_1)

        self.solver = CreateResidualBasedNewtonRaphsonStrategy(
            self.k_model_part, scheme,
            linear_solver, convergence_criteria, builder_and_solver, self.settings["max_iterations"].GetInt(), False,
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

    def IsConverged(self):
        if (hasattr(self, "is_solved_once")):
            return self.GetStrategy().IsConverged()
        else:
            return False

    def SolveCouplingStep(self):
        if (IsBufferInitialized(self)):
            self.is_solved_once = True
            self.solver.Predict()
            self.solver.SolveSolutionStep()
            Kratos.Logger.PrintInfo(self.GetName(), "Solved  formulation.")
            return True

        return False

    def FinializeSolutionStep(self):
        self.solver.FinializeSolutionStep()

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
            elif (scheme_type == "transient"):
                self.is_steady_simulation = False
            else:
                raise Exception("Only \"steady\" and \"transient\" scheme types supported. [ scheme_type = \"" + scheme_type  + "\" ]")
        else:
            raise Exception("\"scheme_type\" is missing in time scheme settings")

    def GetMaxCouplingIterations(self):
        return "N/A"

    def GetModelPart(self):
        return self.k_model_part

    def SetStabilizationMethod(self, stabilization_method):
        if (stabilization_method == "algebraic_flux_corrected"):
            self.element_name = "RansEvmKEpsilonKAFC"
            self.scheme_type = CreateSteadyAlgeraicFluxCorrectedTransportScheme
        elif (stabilization_method == "residual_based_flux_corrected"):
            self.element_name = "RansEvmKEpsilonKResidualBasedFluxCorrected"
            self.scheme_type = CreateSteadyScalarScheme
        elif (stabilization_method == "non_linear_cross_wind_dissipation"):
            self.element_name = "RansEvmKEpsilonKCrossWind"
            self.scheme_type = CreateSteadyScalarScheme
        else:
            raise Exception("Unsupported stabilization method")