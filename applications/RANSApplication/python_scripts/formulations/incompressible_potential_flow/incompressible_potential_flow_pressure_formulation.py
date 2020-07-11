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
from KratosMultiphysics.RANSApplication.formulations.utilities import CreateIncrementalUpdateScheme

class IncompressiblePotentialFlowPressureFormulation(Formulation):
    def __init__(self, model_part, settings):
        super(IncompressiblePotentialFlowPressureFormulation,
              self).__init__(model_part, settings)

        defaults = Kratos.Parameters(r"""{
            "linear_solver_settings": {
                "solver_type": "amgcl"
            },
            "echo_level": 0,
            "relative_tolerance": 1e-12,
            "absolute_tolerance": 1e-12
        }""")

        self.settings.ValidateAndAssignDefaults(defaults)

    def PrepareModelPart(self):
        self.pressure_model_part = CreateFormulationModelPart(
            self, "RansIncompressiblePotentialFlowPressure",
            "RansIncompressiblePotentialFlowPressureBodyForce")

        Kratos.Logger.PrintInfo(self.GetName(),
                                "Created formulation model part.")

    def Initialize(self):
        CalculateNormalsOnConditions(self.GetBaseModelPart())

        solver_settings = self.settings
        linear_solver = CreateLinearSolver(
            solver_settings["linear_solver_settings"])
        builder_and_solver = CreateResidualBasedBlockBuilderAndSolver(
            linear_solver, self.IsPeriodic(), self.GetCommunicator())
        convergence_criteria = CreateResidualCriteria(solver_settings["relative_tolerance"].GetDouble(),
                                                      solver_settings["absolute_tolerance"].GetDouble())
        self.pressure_strategy = CreateResidualBasedNewtonRaphsonStrategy(
            self.pressure_model_part, CreateIncrementalUpdateScheme(),
            linear_solver, convergence_criteria, builder_and_solver, 2, False,
            False, False)

        builder_and_solver.SetEchoLevel(
            solver_settings["echo_level"].GetInt() - 3)
        self.pressure_strategy.SetEchoLevel(
            solver_settings["echo_level"].GetInt() - 2)
        convergence_criteria.SetEchoLevel(
            solver_settings["echo_level"].GetInt() - 1)

        self.pressure_strategy.Initialize()
        Kratos.Logger.PrintInfo(self.GetName(), "Initialized formulation")

    def InitializeSolutionStep(self):
        if (not hasattr(self, "is_initialized")):
            self.is_initialized = True

            RansVariableUtilities.FixFlaggedDofs(self.pressure_model_part,
                                                 KratosRANS.PRESSURE_POTENTIAL,
                                                 Kratos.OUTLET)

            self.pressure_strategy.InitializeSolutionStep()

    def IsConverged(self):
        if (hasattr(self, "is_converged")):
            return self.is_converged
        return False

    def SolveCouplingStep(self):
        if (not self.IsConverged()):
            self.pressure_strategy.Predict()
            self.is_converged = self.pressure_strategy.SolveSolutionStep()
            self.ExecuteAfterCouplingSolveStep()
            Kratos.Logger.PrintInfo(self.GetName(), "Solved  formulation.")
            if (self.is_converged):
                Kratos.Logger.PrintInfo(self.GetName(), "*** CONVERGENCE ACHIEVED ****")
        return True

    def ExecuteAfterCouplingSolveStep(self):
        VariableUtils().CopyModelPartNodalVar(KratosRANS.PRESSURE_POTENTIAL,
                                              Kratos.PRESSURE,
                                              self.pressure_model_part,
                                              self.pressure_model_part, 0)

    def FinializeSolutionStep(self):
        if self.GetBaseModelPart().ProcessInfo[
            Kratos.STEP] + 1 >= 2:
            self.pressure_strategy.FinializeSolutionStep()

    def Check(self):
        self.pressure_strategy.Check()

    def Clear(self):
        self.pressure_strategy.Clear()

    def GetMaxCouplingIterations(self):
        return "N/A"

    def GetModelPart(self):
        return self.pressure_model_part