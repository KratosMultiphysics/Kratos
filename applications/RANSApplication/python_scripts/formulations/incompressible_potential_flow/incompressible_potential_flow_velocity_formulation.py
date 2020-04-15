from __future__ import print_function, absolute_import, division

# import kratos
import KratosMultiphysics as Kratos

# import required applications
import KratosMultiphysics.RANSApplication as KratosRANS

# impot processes
from KratosMultiphysics import IntegrationValuesExtrapolationToNodesProcess as extrapolation_process

# import formulation interface
from KratosMultiphysics.RANSApplication.formulations.formulation import Formulation

# import utilities
from KratosMultiphysics.RANSApplication import RansVariableUtilities
from KratosMultiphysics.RANSApplication.formulations.utilities import CreateLinearSolver
from KratosMultiphysics.RANSApplication.formulations.utilities import CreateFormulationModelPart
from KratosMultiphysics.RANSApplication.formulations.utilities import CalculateNormalsOnConditions
from KratosMultiphysics.RANSApplication.formulations.utilities import CreateResidualBasedBlockBuilderAndSolver
from KratosMultiphysics.RANSApplication.formulations.utilities import CreateResidualCriteria
from KratosMultiphysics.RANSApplication.formulations.utilities import CreateResidualBasedNewtonRaphsonStrategy
from KratosMultiphysics.RANSApplication.formulations.utilities import CreateIncremantalUpdateScheme
from KratosMultiphysics.RANSApplication.formulations.utilities import GetFormulationInfo

class IncompressiblePotentialFlowVelocityFormulation(Formulation):
    def __init__(self, model_part, settings):
        super(IncompressiblePotentialFlowVelocityFormulation, self).__init__(model_part, settings)

        defaults = Kratos.Parameters(r"""{
            "linear_solver_settings": {
                "solver_type": "amgcl"
            },
            "echo_level": 0
        }""")

        self.settings.ValidateAndAssignDefaults(defaults)

    def PrepareModelPart(self):
        self.velocity_model_part = CreateFormulationModelPart(self,
                                                                "RansIncompressiblePotentialFlowVelocity",
                                                                "RansIncompressiblePotentialFlowVelocity")

        Kratos.Logger.PrintInfo(self.GetName(), "Created formulation model part.")

    def Initialize(self):
        CalculateNormalsOnConditions(self.GetBaseModelPart())

        solver_settings = self.settings
        linear_solver = CreateLinearSolver(solver_settings["linear_solver_settings"])
        builder_and_solver = CreateResidualBasedBlockBuilderAndSolver(linear_solver, self.IsPeriodic(), self.GetCommunicator())
        convergence_criteria = CreateResidualCriteria(1e-12, 1e-12)
        self.velocity_strategy = CreateResidualBasedNewtonRaphsonStrategy(
                    self.velocity_model_part, CreateIncremantalUpdateScheme(), linear_solver,
                    convergence_criteria, builder_and_solver, 2, False, False, False)

        builder_and_solver.SetEchoLevel(
            solver_settings["echo_level"].GetInt() - 3)
        self.velocity_strategy.SetEchoLevel(
            solver_settings["echo_level"].GetInt() - 2)
        convergence_criteria.SetEchoLevel(
            solver_settings["echo_level"].GetInt() - 1)

        self.velocity_strategy.Initialize()
        Kratos.Logger.PrintInfo(self.GetName(), "Initialized formulation")

    def InitializeSolutionStep(self):
        if (not hasattr(self, "is_initialized")):
            self.is_initialized = True

            RansVariableUtilities.FixFlaggedDofs(self.velocity_model_part,
                                                 KratosRANS.VELOCITY_POTENTIAL,
                                                 Kratos.OUTLET)

            self.velocity_strategy.InitializeSolutionStep()

    def IsConverged(self):
        if (hasattr(self, "is_solved")):
            return self.is_solved
        return False

    def SolveCouplingStep(self):
        self.is_solved = True
        self.velocity_strategy.Predict()
        self.velocity_strategy.SolveSolutionStep()

    def ExecuteAfterCouplingSolveStep(self):
        RansVariableUtilities.CopyFlaggedVariableToNonHistorical(
            self.velocity_model_part, Kratos.VELOCITY, Kratos.INLET)

        # extrapolate gauss point velocities to nodal velocities
        extrapolation_settings = Kratos.Parameters('''
        {
            "echo_level"                 : 0,
            "area_average"               : true,
            "average_variable"           : "NODAL_AREA",
            "list_of_variables"          : ["VELOCITY"],
            "extrapolate_non_historical" : false
        }''')
        extrapolation_process(self.velocity_model_part,
                                extrapolation_settings).Execute()

        # take back the original inlet velocities
        RansVariableUtilities.CopyFlaggedVariableFromNonHistorical(
            self.velocity_model_part, Kratos.VELOCITY, Kratos.INLET)

    def FinializeSolutionStep(self):
        self.velocity_strategy.FinializeSolutionStep()

    def Check(self):
        self.velocity_strategy.Check()

    def Clear(self):
        self.velocity_strategy.Clear()

    def GetInfo(self):
        return GetFormulationInfo(self, self.velocity_model_part)

