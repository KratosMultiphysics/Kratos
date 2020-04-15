from __future__ import print_function, absolute_import, division

# import kratos
import KratosMultiphysics as Kratos

# import required applications
import KratosMultiphysics.RANSApplication as KratosRANS

# import formulation interface
from KratosMultiphysics.RANSApplication.formulations.formulation import Formulation

# import utilities
from KratosMultiphysics import VariableUtils
from KratosMultiphysics.RANSApplication import RansVariableUtilities
from KratosMultiphysics.RANSApplication.formulations.utilities import CreateLinearSolver
from KratosMultiphysics.RANSApplication.formulations.utilities import CreateFormulationModelPart
from KratosMultiphysics.RANSApplication.formulations.utilities import CalculateNormalsOnConditions
from KratosMultiphysics.RANSApplication.formulations.utilities import CreateResidualBasedBlockBuilderAndSolver
from KratosMultiphysics.RANSApplication.formulations.utilities import CreateResidualCriteria
from KratosMultiphysics.RANSApplication.formulations.utilities import CreateResidualBasedNewtonRaphsonStrategy
from KratosMultiphysics.RANSApplication.formulations.utilities import CreateSteadyScalarScheme
from KratosMultiphysics.RANSApplication.formulations.utilities import CreateBossakScalarScheme
from KratosMultiphysics.RANSApplication.formulations.utilities import GetFormulationInfo

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
            },
            "reform_dofs_at_each_step": true,
            "move_mesh_strategy": 0,
            "move_mesh_flag": false,
            "compute_reactions": false
        }""")

        self.settings.ValidateAndAssignDefaults(defaults)
        self.echo_level = self.settings["echo_level"].GetInt()

    def PrepareModelPart(self):
        if self.GetBaseModelPart().ProcessInfo[Kratos.DOMAIN_SIZE] == 2:
            condition_name = "LineCondition"
        else:
            condition_name = "SurfaceCondition"

        self.k_model_part = CreateFormulationModelPart(self, "RansEvmKEpsilonK", condition_name)

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
            scheme = CreateSteadyScalarScheme(self.settings["relaxation_factor"].GetDouble())
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
        if (hasattr(self, "is_solved")):
            return self.is_solved
        return False

    def SolveCouplingStep(self):
        self.solver.Predict()
        self.is_solved = self.solver.SolveSolutionStep()

    def FinializeSolutionStep(self):
        self.solver.FinializeSolutionStep()

    def Check(self):
        self.solver.Check()

    def Clear(self):
        self.solver.Clear()

    def GetInfo(self):
        return GetFormulationInfo(self, self.k_model_part)

    def GetStrategy(self):
        return self.solver