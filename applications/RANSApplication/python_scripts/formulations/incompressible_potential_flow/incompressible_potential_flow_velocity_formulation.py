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
from KratosMultiphysics.RANSApplication.formulations.utilities import CreateIncrementalUpdateScheme

class IncompressiblePotentialFlowVelocityFormulation(Formulation):
    def __init__(self, model_part, settings):
        super(IncompressiblePotentialFlowVelocityFormulation, self).__init__(model_part, settings)

        defaults = Kratos.Parameters(r"""{
            "linear_solver_settings": {
                "solver_type": "amgcl"
            },
            "echo_level": 0,
            "relative_tolerance": 1e-12,
            "absolute_tolerance": 1e-12
        }""")

        self.settings.ValidateAndAssignDefaults(defaults)
        self.max_coupling_iterations = 1

    def PrepareModelPart(self):
        self.velocity_model_part = CreateFormulationModelPart(self,
                                                                "RansIncompressiblePotentialFlowVelocity",
                                                                "RansIncompressiblePotentialFlowVelocityInlet")

        Kratos.Logger.PrintInfo(self.GetName(), "Created formulation model part.")

    def Initialize(self):
        CalculateNormalsOnConditions(self.GetBaseModelPart())

        solver_settings = self.settings
        linear_solver = CreateLinearSolver(solver_settings["linear_solver_settings"])
        builder_and_solver = CreateResidualBasedBlockBuilderAndSolver(linear_solver, self.IsPeriodic(), self.GetCommunicator())
        convergence_criteria = CreateResidualCriteria(
            [(KratosRANS.VELOCITY_POTENTIAL, solver_settings["relative_tolerance"].GetDouble(
            ), solver_settings["absolute_tolerance"].GetDouble())])
        self.velocity_strategy = CreateResidualBasedNewtonRaphsonStrategy(
                    self.velocity_model_part, CreateIncrementalUpdateScheme(),
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

            Kratos.VariableUtils().ApplyFixity(KratosRANS.VELOCITY_POTENTIAL,
                                               True, self.velocity_model_part.Nodes, Kratos.OUTLET, True)

            self.velocity_strategy.InitializeSolutionStep()

    def IsConverged(self):
        if (hasattr(self, "is_converged")):
            return self.is_converged
        return False

    def SolveCouplingStep(self):
        if (not self.IsConverged()):
            self.velocity_strategy.Predict()
            self.is_converged = self.velocity_strategy.SolveSolutionStep()
            self.ExecuteAfterCouplingSolveStep()
            Kratos.Logger.PrintInfo(self.GetName(), "Solved  formulation.")
            if (self.is_converged):
                Kratos.Logger.PrintInfo(self.GetName(), "*** CONVERGENCE ACHIEVED ****")
        return True

    def ExecuteAfterCouplingSolveStep(self):
        Kratos.VariableUtils().CopyModelPartFlaggedNodalHistoricalVarToNonHistoricalVar(
            Kratos.VELOCITY, self.velocity_model_part, self.velocity_model_part, Kratos.INLET, True, 0)

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
        Kratos.VariableUtils().CopyModelPartFlaggedNodalNonHistoricalVarToHistoricalVar(
            Kratos.VELOCITY, self.velocity_model_part, self.velocity_model_part, Kratos.INLET, True, 0)

    def FinializeSolutionStep(self):
        self.velocity_strategy.FinializeSolutionStep()

    def Check(self):
        self.velocity_strategy.Check()

    def Clear(self):
        self.velocity_strategy.Clear()

    def GetMaxCouplingIterations(self):
        return "N/A"

    def GetModelPart(self):
        return self.velocity_model_part
