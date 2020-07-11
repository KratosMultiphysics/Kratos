# Importing the Kratos Library
import KratosMultiphysics as Kratos
from KratosMultiphysics import IsDistributedRun
from KratosMultiphysics.kratos_utilities import CheckIfApplicationsAvailable

# Import applications
import KratosMultiphysics.FluidDynamicsApplication as KratosCFD
import KratosMultiphysics.RANSApplication as KratosRANS

# import formulation interface
from KratosMultiphysics.RANSApplication.formulations.formulation import Formulation

# import utilities
from KratosMultiphysics.RANSApplication import RansCalculationUtilities
from KratosMultiphysics.RANSApplication.formulations.utilities import CreateLinearSolver
from KratosMultiphysics.RANSApplication.formulations.utilities import CreateResidualBasedNewtonRaphsonStrategy
from KratosMultiphysics.RANSApplication.formulations.utilities import CreateResidualBasedBlockBuilderAndSolver
from KratosMultiphysics.RANSApplication.formulations.utilities import CreateResidualCriteria
from KratosMultiphysics.RANSApplication.formulations.utilities import CreateFormulationModelPart
from KratosMultiphysics.RANSApplication.formulations.utilities import CalculateNormalsOnConditions
from KratosMultiphysics.RANSApplication.formulations.utilities import IsBufferInitialized
from KratosMultiphysics.RANSApplication.formulations.utilities import InitializePeriodicConditions

# case specific imports
if (IsDistributedRun() and CheckIfApplicationsAvailable("TrilinosApplication")):
    from KratosMultiphysics.TrilinosApplication import TrilinosUPCriteria as up_criteria
    from KratosMultiphysics.FluidDynamicsApplication.TrilinosExtension import TrilinosPredictorCorrectorVelocityBossakSchemeTurbulent as dynamic_scheme
    from KratosMultiphysics.FluidDynamicsApplication.TrilinosExtension import TrilinosResidualBasedSimpleSteadyScheme as steady_scheme
elif (not IsDistributedRun()):
    from KratosMultiphysics.FluidDynamicsApplication import VelPrCriteria as up_criteria
    from KratosMultiphysics.FluidDynamicsApplication import ResidualBasedPredictorCorrectorVelocityBossakSchemeTurbulent as dynamic_scheme
    from KratosMultiphysics.FluidDynamicsApplication import ResidualBasedSimpleSteadyScheme as steady_scheme
else:
    raise Exception("Distributed run requires TrilinosApplication")


class MonolithicVelocityPressureFormulation(Formulation):
    def __init__(self, model_part, settings):
        super(MonolithicVelocityPressureFormulation, self).__init__(model_part, settings)

        ##settings string in json format
        default_settings = Kratos.Parameters("""
        {
            "formulation_name": "monolithic",
            "maximum_iterations": 10,
            "echo_level": 0,
            "compute_reactions": false,
            "reform_dofs_at_each_step": false,
            "relative_velocity_tolerance": 1e-5,
            "absolute_velocity_tolerance": 1e-7,
            "relative_pressure_tolerance": 1e-5,
            "absolute_pressure_tolerance": 1e-7,
            "linear_solver_settings": {
                "solver_type": "amgcl"
            },
            "move_mesh_strategy": 0,
            "move_mesh_flag": false,
            "velocity_relaxation":0.9,
            "pressure_relaxation":0.9,
            "oss_switch": 0,
            "dynamic_tau": 0.01
        }""")

        settings.ValidateAndAssignDefaults(default_settings)

        self.min_buffer_size = 2
        self.element_has_nodal_properties = True
        self.fractional_step_model_part = None

        ## Construct the linear solvers
        self.linear_solver = CreateLinearSolver(self.settings["linear_solver_settings"])
        self.echo_level = self.settings["echo_level"].GetInt()

        self.compute_reactions = self.settings["compute_reactions"].GetBool()
        self.SetMaxCouplingIterations(1)

        Kratos.Logger.PrintInfo(self.GetName(), "Construction of formulation finished.")

    def AddVariables(self):
        self.GetBaseModelPart().AddNodalSolutionStepVariable(Kratos.VELOCITY)
        self.GetBaseModelPart().AddNodalSolutionStepVariable(Kratos.ACCELERATION)
        self.GetBaseModelPart().AddNodalSolutionStepVariable(Kratos.MESH_VELOCITY)
        self.GetBaseModelPart().AddNodalSolutionStepVariable(Kratos.PRESSURE)
        self.GetBaseModelPart().AddNodalSolutionStepVariable(Kratos.IS_STRUCTURE)
        self.GetBaseModelPart().AddNodalSolutionStepVariable(Kratos.DISPLACEMENT)
        self.GetBaseModelPart().AddNodalSolutionStepVariable(Kratos.VISCOSITY)
        self.GetBaseModelPart().AddNodalSolutionStepVariable(Kratos.DENSITY)
        self.GetBaseModelPart().AddNodalSolutionStepVariable(Kratos.BODY_FORCE)
        self.GetBaseModelPart().AddNodalSolutionStepVariable(Kratos.NODAL_AREA)
        self.GetBaseModelPart().AddNodalSolutionStepVariable(Kratos.NODAL_H)
        self.GetBaseModelPart().AddNodalSolutionStepVariable(Kratos.ADVPROJ)
        self.GetBaseModelPart().AddNodalSolutionStepVariable(Kratos.DIVPROJ)
        self.GetBaseModelPart().AddNodalSolutionStepVariable(Kratos.REACTION)
        self.GetBaseModelPart().AddNodalSolutionStepVariable(Kratos.REACTION_WATER_PRESSURE)
        self.GetBaseModelPart().AddNodalSolutionStepVariable(Kratos.EXTERNAL_PRESSURE)
        self.GetBaseModelPart().AddNodalSolutionStepVariable(Kratos.NORMAL)
        self.GetBaseModelPart().AddNodalSolutionStepVariable(Kratos.Y_WALL)
        self.GetBaseModelPart().AddNodalSolutionStepVariable(KratosCFD.Q_VALUE)
        self.GetBaseModelPart().AddNodalSolutionStepVariable(Kratos.KINEMATIC_VISCOSITY)
        self.GetBaseModelPart().AddNodalSolutionStepVariable(KratosRANS.TURBULENT_KINETIC_ENERGY)

        Kratos.Logger.PrintInfo(self.GetName(), "Added solution step variables.")

    def AddDofs(self):
        Kratos.VariableUtils().AddDof(Kratos.VELOCITY_X, Kratos.REACTION_X,self.GetBaseModelPart())
        Kratos.VariableUtils().AddDof(Kratos.VELOCITY_Y, Kratos.REACTION_Y,self.GetBaseModelPart())
        Kratos.VariableUtils().AddDof(Kratos.VELOCITY_Z, Kratos.REACTION_Z,self.GetBaseModelPart())
        Kratos.VariableUtils().AddDof(Kratos.PRESSURE, Kratos.REACTION_WATER_PRESSURE,self.GetBaseModelPart())

        Kratos.Logger.PrintInfo(self.GetName(), "Added dofs.")

    def PrepareModelPart(self):
        self.monolithic_model_part = CreateFormulationModelPart(self,
                                                                "VMS",
                                                                self.condition_name)
        Kratos.Logger.PrintInfo(self.GetName(), "Created formulation model part.")

    def Initialize(self):
        model_part = self.GetBaseModelPart()
        CalculateNormalsOnConditions(model_part)

        process_info = model_part.ProcessInfo
        bossak_alpha = process_info[Kratos.BOSSAK_ALPHA]

        if (self.IsPeriodic()):
            if (self.domain_size == 2):
                periodic_variables_list = [Kratos.VELOCITY_X, Kratos.VELOCITY_Y, Kratos.PRESSURE]
            else:
                periodic_variables_list = [Kratos.VELOCITY_X, Kratos.VELOCITY_Y, Kratos.VELOCITY_Z, Kratos.PRESSURE]
            InitializePeriodicConditions(model_part,
                                         self.monolithic_model_part,
                                         periodic_variables_list)

        if self.is_steady_simulation:
            conv_criteria = CreateResidualCriteria(self.settings["relative_velocity_tolerance"].GetDouble(),
                                                   self.settings["absolute_velocity_tolerance"].GetDouble())
            scheme = steady_scheme(self.settings["velocity_relaxation"].GetDouble(),
                                   self.settings["pressure_relaxation"].GetDouble(),
                                   self.domain_size)
        else:
            conv_criteria = up_criteria(self.settings["relative_velocity_tolerance"].GetDouble(),
                                        self.settings["absolute_velocity_tolerance"].GetDouble(),
                                        self.settings["relative_pressure_tolerance"].GetDouble(),
                                        self.settings["absolute_pressure_tolerance"].GetDouble())
            scheme = dynamic_scheme(bossak_alpha,
                                    self.settings["move_mesh_strategy"].GetInt(),
                                    self.domain_size)

        builder_and_solver = CreateResidualBasedBlockBuilderAndSolver(
                                    self.linear_solver,
                                    self.IsPeriodic(),
                                    self.GetCommunicator())

        self.solver = CreateResidualBasedNewtonRaphsonStrategy(self.monolithic_model_part,
                                                               scheme,
                                                               self.linear_solver,
                                                               conv_criteria,
                                                               builder_and_solver,
                                                               self.settings["maximum_iterations"].GetInt(),
                                                               self.settings["compute_reactions"].GetBool(),
                                                               self.settings["reform_dofs_at_each_step"].GetBool(),
                                                               self.settings["move_mesh_flag"].GetBool())

        builder_and_solver.SetEchoLevel(
            self.settings["echo_level"].GetInt() - 3)
        self.solver.SetEchoLevel(
            self.settings["echo_level"].GetInt() - 2)
        conv_criteria.SetEchoLevel(
            self.settings["echo_level"].GetInt() - 1)

        self.GetBaseModelPart().ProcessInfo.SetValue(Kratos.DYNAMIC_TAU, self.settings["dynamic_tau"].GetDouble())
        self.GetBaseModelPart().ProcessInfo.SetValue(Kratos.OSS_SWITCH, self.settings["oss_switch"].GetInt())

        super(MonolithicVelocityPressureFormulation, self).Initialize()
        self.solver.Initialize()

        Kratos.Logger.PrintInfo(self.GetName(), "Solver initialization finished.")

    def GetMinimumBufferSize(self):
        return self.min_buffer_size

    def Finalize(self):
        self.solver.Clear()
        super(MonolithicVelocityPressureFormulation, self).Finalize()

    def SolveCouplingStep(self):
        if (IsBufferInitialized(self)):
            max_iterations = self.GetMaxCouplingIterations()
            for iteration in range(max_iterations):
                self.solver.Predict()
                self.is_converged = self.solver.SolveSolutionStep()
                self.ExecuteAfterCouplingSolveStep()
                Kratos.Logger.PrintInfo(self.GetName(), "Solved coupling iteration " + str(iteration + 1) + "/" + str(max_iterations) + ".")
                return True

        return False

    def InitializeSolutionStep(self):
        if (IsBufferInitialized(self)):
            super(MonolithicVelocityPressureFormulation, self).InitializeSolutionStep()
            self.solver.InitializeSolutionStep()

    def FinializeSolutionStep(self):
        if (IsBufferInitialized(self)):
            self.solver.FinializeSolutionStep()
            super(MonolithicVelocityPressureFormulation, self).FinializeSolutionStep()

    def Check(self):
        super(MonolithicVelocityPressureFormulation, self).Check()
        self.solver.Check()

    def Clear(self):
        self.solver.Clear()
        super(MonolithicVelocityPressureFormulation, self).Clear()

    def SetTimeSchemeSettings(self, settings):
        if (settings.Has("scheme_type")):
            scheme_type = settings["scheme_type"].GetString()
            if (scheme_type == "steady"):
                self.is_steady_simulation = True
            elif (scheme_type == "transient"):
                self.is_steady_simulation = False
                default_settings = Kratos.Parameters('''{
                    "scheme_type": "transient",
                    "alpha_bossak": -0.3
                }''')
                settings.ValidateAndAssignDefaults(default_settings)
                self.GetBaseModelPart().ProcessInfo.SetValue(Kratos.BOSSAK_ALPHA, settings["alpha_bossak"].GetDouble())
            else:
                raise Exception("Only \"steady\" and \"transient\" scheme types supported. [ scheme_type = \"" + scheme_type  + "\" ]")
        else:
            raise Exception("\"scheme_type\" is missing in time scheme settings")

        self.time_scheme_settings = settings

    def SetConstants(self, settings):
        defaults = Kratos.Parameters('''{
            "von_karman": 0.41,
            "beta"      : 5.2,
            "c_mu"      : 0.09
        }''')
        settings.ValidateAndAssignDefaults(defaults)

        # set constants
        self.von_karman = settings["von_karman"].GetDouble()
        self.beta = settings["beta"].GetDouble()
        self.y_plus_limit = RansCalculationUtilities.CalculateLogarithmicYPlusLimit(
                                                                                self.von_karman,
                                                                                self.beta
                                                                                )

        self.GetBaseModelPart().ProcessInfo.SetValue(KratosRANS.WALL_VON_KARMAN, self.von_karman)
        self.GetBaseModelPart().ProcessInfo.SetValue(KratosRANS.WALL_SMOOTHNESS_BETA, self.beta)
        self.GetBaseModelPart().ProcessInfo.SetValue(KratosRANS.RANS_LINEAR_LOG_LAW_Y_PLUS_LIMIT, self.y_plus_limit)
        self.GetBaseModelPart().ProcessInfo.SetValue(KratosRANS.TURBULENCE_RANS_C_MU, settings["c_mu"].GetDouble())

    def SetWallFunctionSettings(self, settings):
        wall_function_region_type = "logarithmic_region_only"
        if (settings.Has("wall_function_region_type")):
            wall_function_region_type = settings[
                "wall_function_region_type"].GetString()

        if (wall_function_region_type == "logarithmic_region_only"):
            self.condition_name = "RansVMSMonolithicKBasedWall"
        else:
            msg = "Unsupported wall function region type provided. [ wall_function_region_type = \"" + wall_function_region_type + "\" ]."
            msg += "Supported wall function region types are:\n"
            msg += "\tlogarithmic_region_only\n"
            msg += "\tlinear_and_logarithmic_regions\n"
            raise Exception(msg)

    def GetStrategy(self):
        return self.solver

    def GetModelPart(self):
        return self.monolithic_model_part
