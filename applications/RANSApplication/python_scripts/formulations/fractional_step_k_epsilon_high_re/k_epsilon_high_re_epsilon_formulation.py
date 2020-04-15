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

class KEpsilonHighReEpsilonFormulation(Formulation):
    def __init__(self, model_part, settings):
        super(KEpsilonHighReEpsilonFormulation,
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
        self.epsilon_model_part = CreateFormulationModelPart(self, "RansEvmKEpsilonEpsilon", "RansEvmKEpsilonEpsilonWall")

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
                self.epsilon_model_part.ProcessInfo[Kratos.BOSSAK_ALPHA],
                self.settings["relaxation_factor"].GetDouble(),
                KratosRANS.TURBULENT_ENERGY_DISSIPATION_RATE,
                KratosRANS.TURBULENT_ENERGY_DISSIPATION_RATE_2,
                KratosRANS.RANS_AUXILIARY_VARIABLE_2)

        self.solver = CreateResidualBasedNewtonRaphsonStrategy(
            self.epsilon_model_part, scheme,
            linear_solver, convergence_criteria, builder_and_solver, self.settings["max_iterations"].GetInt(), False,
            False, False)

        builder_and_solver.SetEchoLevel(
            solver_settings["echo_level"].GetInt() - 3)
        self.solver.SetEchoLevel(
            solver_settings["echo_level"].GetInt() - 2)
        convergence_criteria.SetEchoLevel(
            solver_settings["echo_level"].GetInt() - 1)

        process_info = self.epsilon_model_part.ProcessInfo

        settings = Kratos.Parameters("{" + """
            "model_part_name" : "{0:s}",
            "echo_level"      : {1:d},
            "c_mu"            : {2:f},
            "min_value"       : 1e-15
        """.format(self.epsilon_model_part.Name, self.echo_level-1, process_info[KratosRANS.TURBULENCE_RANS_C_MU]) + "}")
        self.nut_process = KratosRANS.RansNutKEpsilonHighReCalculationProcess(self.epsilon_model_part.GetModel(), settings)

        settings = Kratos.Parameters("{" + """
            "model_part_name" : "{0:s}",
            "echo_level"      : {1:d},
            "c_mu"            : {2:f},
            "von_karman"      : {3:f},
            "beta"            : {4:f},
            "min_value"       : 1e-15
        """.format(self.epsilon_model_part.Name,
                    self.echo_level-1,
                    process_info[KratosRANS.TURBULENCE_RANS_C_MU],
                    process_info[KratosRANS.WALL_VON_KARMAN],
                    process_info[KratosRANS.WALL_SMOOTHNESS_BETA]) + "}")
        self.nut_wall_process = KratosRANS.RansNutYPlusWallFunctionProcess(self.epsilon_model_part.GetModel(), settings)

        self.nut_process.ExecuteInitialize()
        self.nut_wall_process.ExecuteInitialize()

        self.solver.Initialize()
        Kratos.Logger.PrintInfo(self.GetName(), "Initialized formulation")

    def InitializeSolutionStep(self):
        self.nut_process.ExecuteInitializeSolutionStep()
        self.nut_wall_process.ExecuteInitializeSolutionStep()
        self.solver.InitializeSolutionStep()

    def IsConverged(self):
        if (hasattr(self, "is_solved")):
            return self.is_solved
        return False

    def SolveCouplingStep(self):
        self.solver.Predict()
        self.is_solved = self.solver.SolveSolutionStep()

    def ExecuteAfterCouplingSolveStep(self):
        self.nut_process.Execute()
        self.nut_wall_process.Execute()

    def FinializeSolutionStep(self):
        self.solver.FinializeSolutionStep()
        self.nut_process.ExecuteFinializeSolutionStep()
        self.nut_wall_process.ExecuteFinializeSolutionStep()

    def Check(self):
        self.nut_process.Check()
        self.nut_wall_process.Check()
        self.solver.Check()

    def Clear(self):
        self.solver.Clear()

    def GetInfo(self):
        return GetFormulationInfo(self, self.epsilon_model_part)

    def GetStrategy(self):
        return self.solver

