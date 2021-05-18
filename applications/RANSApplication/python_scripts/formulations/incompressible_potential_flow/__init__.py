# import kratos
import KratosMultiphysics as Kratos

# import RANS
import KratosMultiphysics.RANSApplication as KratosRANS

# impot processes
from KratosMultiphysics import IntegrationValuesExtrapolationToNodesProcess as extrapolation_process

# import formulation interface
from KratosMultiphysics.RANSApplication.formulations.rans_formulation import RansFormulation

from KratosMultiphysics.RANSApplication.formulations.utilities import CreateRansFormulationModelPart
from KratosMultiphysics.RANSApplication.formulations.utilities import CalculateNormalsOnConditions
from KratosMultiphysics.RANSApplication.formulations.utilities import CreateBlockBuilderAndSolver
from KratosMultiphysics.RANSApplication.formulations.utilities import GetKratosObjectPrototype

class IncompressiblePotentialFlowRansFormulation(RansFormulation):
    def __init__(self, model_part, settings):
        """Incompressible potential flow rans formulation

        This RansFormulation solves incompressible potential flow equation for steady
        state problems. (The dof is KratosRANS.VELOCITY_POTENTIAL) Afterwards, Kratos.VELOCITY
        variable is updated properly from computed KratosRANS.VELOCITY_POTENTIAL values for each node.

        This does not support periodic conditions.

        Args:
            model_part (Kratos.ModelPart): ModelPart to be used in the formulation.
            settings (Kratos.Parameters): Settings to be used in the formulation.
        """
        super().__init__(model_part, settings)

        default_settings = Kratos.Parameters(r'''
        {
            "formulation_name": "incompressible_potential_flow",
            "linear_solver_settings": {
                "solver_type": "amgcl"
            },
            "echo_level": 0,
            "relative_tolerance": 1e-12,
            "absolute_tolerance": 1e-12
        }''')
        self.GetParameters().ValidateAndAssignDefaults(default_settings)
        self.SetMaxCouplingIterations(1)

    def AddVariables(self):
        self.GetBaseModelPart().AddNodalSolutionStepVariable(Kratos.VELOCITY)
        self.GetBaseModelPart().AddNodalSolutionStepVariable(Kratos.PRESSURE)
        self.GetBaseModelPart().AddNodalSolutionStepVariable(KratosRANS.VELOCITY_POTENTIAL)
        self.GetBaseModelPart().AddNodalSolutionStepVariable(Kratos.BODY_FORCE)
        self.GetBaseModelPart().AddNodalSolutionStepVariable(Kratos.DENSITY)
        self.GetBaseModelPart().AddNodalSolutionStepVariable(Kratos.NORMAL)

        # these variables are required because of outlet boundary condition
        self.GetBaseModelPart().AddNodalSolutionStepVariable(Kratos.EXTERNAL_PRESSURE)

        # these variables are required for calculation of reactions
        self.GetBaseModelPart().AddNodalSolutionStepVariable(Kratos.REACTION)
        self.GetBaseModelPart().AddNodalSolutionStepVariable(Kratos.REACTION_WATER_PRESSURE)

        Kratos.Logger.PrintInfo(self.__class__.__name__, "Added solution step variables.")

    def AddDofs(self):
        Kratos.VariableUtils().AddDof(KratosRANS.VELOCITY_POTENTIAL, self.GetBaseModelPart())

        # these dofs are required since they are fixed in inlet and outlet boundary conditions.
        Kratos.VariableUtils().AddDof(Kratos.VELOCITY_X, Kratos.REACTION_X,self.GetBaseModelPart())
        Kratos.VariableUtils().AddDof(Kratos.VELOCITY_Y, Kratos.REACTION_Y,self.GetBaseModelPart())
        Kratos.VariableUtils().AddDof(Kratos.VELOCITY_Z, Kratos.REACTION_Z,self.GetBaseModelPart())
        Kratos.VariableUtils().AddDof(Kratos.PRESSURE, Kratos.REACTION_WATER_PRESSURE,self.GetBaseModelPart())

        Kratos.Logger.PrintInfo(self.__class__.__name__, "Added solution step dofs.")

    def GetMinimumBufferSize(self):
        return 1

    def IsPeriodic(self):
        return False

    def SetIsPeriodic(self, value):
        if (value):
            raise Exception("Periodic conditions are not supported by incompressible potential flow solver.")

    def PrepareModelPart(self):
        self.velocity_model_part = CreateRansFormulationModelPart(
            self.GetComputingModelPart(),
            self.__class__.__name__,
            self.GetDomainSize(),
            "RansIncompressiblePotentialFlowVelocity",
            "RansIncompressiblePotentialFlowVelocityInlet")

        Kratos.Logger.PrintInfo(self.__class__.__name__, "Created formulation model part.")

    def Initialize(self):
        CalculateNormalsOnConditions(self.GetBaseModelPart())

        solver_settings = self.GetParameters()

        linear_solver_factory = GetKratosObjectPrototype("LinearSolverFactory")
        linear_solver = linear_solver_factory(solver_settings["linear_solver_settings"])

        builder_and_solver = CreateBlockBuilderAndSolver(
            linear_solver,
            self.IsPeriodic(),
            self.GetCommunicator())

        convergence_criteria_type = GetKratosObjectPrototype("MixedGenericCriteria")
        convergence_criteria = convergence_criteria_type(
            [(KratosRANS.VELOCITY_POTENTIAL, solver_settings["relative_tolerance"].GetDouble(),
            solver_settings["absolute_tolerance"].GetDouble())])

        scheme_type = GetKratosObjectPrototype("ResidualBasedIncrementalUpdateStaticScheme")
        strategy_type = GetKratosObjectPrototype("ResidualBasedNewtonRaphsonStrategy")

        self.velocity_strategy = strategy_type(
            self.velocity_model_part, scheme_type(),
            convergence_criteria, builder_and_solver, 2, False, False, False)

        builder_and_solver.SetEchoLevel(solver_settings["echo_level"].GetInt())
        self.velocity_strategy.SetEchoLevel(solver_settings["echo_level"].GetInt())
        convergence_criteria.SetEchoLevel(solver_settings["echo_level"].GetInt())

        super().Initialize()
        Kratos.Logger.PrintInfo(self.__class__.__name__, "Initialized formulation")

    def InitializeSolutionStep(self):
        if (not hasattr(self, "is_initialized")):
            self.is_initialized = True

            Kratos.VariableUtils().ApplyFixity(
                KratosRANS.VELOCITY_POTENTIAL,
                True,
                self.velocity_model_part.Nodes,
                Kratos.OUTLET,
                True)

            super().InitializeSolutionStep()

    def IsConverged(self):
        if (hasattr(self, "is_converged")):
            return self.is_converged
        return False

    def SolveCouplingStep(self):
        if (not self.IsConverged()):
            self.velocity_strategy.Predict()
            self.is_converged = self.velocity_strategy.SolveSolutionStep()
            self.ExecuteAfterCouplingSolveStep()
            Kratos.Logger.PrintInfo(self.__class__.__name__, "Solved  formulation.")
            if (self.is_converged):
                Kratos.Logger.PrintInfo(self.__class__.__name__, "*** CONVERGENCE ACHIEVED ****")
        return True

    def ExecuteAfterCouplingSolveStep(self):
        Kratos.VariableUtils().CopyModelPartFlaggedNodalHistoricalVarToNonHistoricalVar(
            Kratos.VELOCITY,
            self.velocity_model_part,
            self.velocity_model_part,
            Kratos.INLET,
            True,
            0)

        # extrapolate gauss point velocities to nodal velocities
        extrapolation_settings = Kratos.Parameters('''
        {
            "echo_level"                 : 0,
            "area_average"               : true,
            "average_variable"           : "NODAL_AREA",
            "list_of_variables"          : ["VELOCITY"],
            "extrapolate_non_historical" : false
        }''')
        extrapolation_process(
            self.velocity_model_part,
            extrapolation_settings).Execute()

        # take back the original inlet velocities
        Kratos.VariableUtils().CopyModelPartFlaggedNodalNonHistoricalVarToHistoricalVar(
            Kratos.VELOCITY,
            self.velocity_model_part,
            self.velocity_model_part,
            Kratos.INLET,
            True,
            0)

    def GetStrategy(self):
        return self.velocity_strategy

    def GetMaxCouplingIterations(self):
        return 0
