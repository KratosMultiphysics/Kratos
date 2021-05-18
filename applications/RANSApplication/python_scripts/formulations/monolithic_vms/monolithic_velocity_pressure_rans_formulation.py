# Importing the Kratos Library
import KratosMultiphysics as Kratos

# Import applications
import KratosMultiphysics.FluidDynamicsApplication as KratosCFD
import KratosMultiphysics.RANSApplication as KratosRANS

# import formulation interface
from KratosMultiphysics.RANSApplication.formulations.rans_formulation import RansFormulation

# import utilities
from KratosMultiphysics.RANSApplication.formulations.utilities import CreateBlockBuilderAndSolver
from KratosMultiphysics.RANSApplication.formulations.utilities import CreateRansFormulationModelPart
from KratosMultiphysics.RANSApplication.formulations.utilities import CalculateNormalsOnConditions
from KratosMultiphysics.RANSApplication.formulations.utilities import InitializePeriodicConditions
from KratosMultiphysics.RANSApplication.formulations.utilities import GetKratosObjectPrototype

class StabilizedFormulation(object):
    """Helper class to define stabilization-dependent parameters."""
    def __init__(self, settings):
        self.element_name = None
        self.element_has_nodal_properties = False
        self.process_data = {}

        if settings.Has("element_type"):
            element_type = settings["element_type"].GetString()
            if element_type == "vms":
                self._SetUpClassicVMS(settings)
            elif element_type == "qsvms":
                self._SetUpQSVMS(settings)
            else:
                raise RuntimeError("Unknown element_type. [ element_type = \"" + element_type + "\"")
        else:
            print(settings)
            raise RuntimeError("Argument \'element_type\' not found in stabilization settings.")

    def SetProcessInfo(self, model_part):
        for variable,value in self.process_data.items():
            model_part.ProcessInfo[variable] = value

    def _SetUpClassicVMS(self, settings):
        default_settings = Kratos.Parameters(r"""{
            "element_type": "vms",
            "use_orthogonal_subscales": false,
            "dynamic_tau": 0.01
        }""")

        self.element_name = 'VMS'
        settings.ValidateAndAssignDefaults(default_settings)

        # set the nodal material properties flag
        self.element_has_nodal_properties = True

        self.process_data[Kratos.DYNAMIC_TAU] = settings["dynamic_tau"].GetDouble()
        use_oss = settings["use_orthogonal_subscales"].GetBool()
        self.process_data[Kratos.OSS_SWITCH] = int(use_oss)

    def _SetUpQSVMS(self, settings):
        default_settings = Kratos.Parameters(r"""{
            "element_type": "qsvms",
            "use_orthogonal_subscales": false,
            "dynamic_tau": 0.0,
            "element_manages_time_integration": false
        }""")
        settings.ValidateAndAssignDefaults(default_settings)

        self.element_name = "QSVMS"

        self.process_data[Kratos.DYNAMIC_TAU] = settings["dynamic_tau"].GetDouble()
        use_oss = settings["use_orthogonal_subscales"].GetBool()
        self.process_data[Kratos.OSS_SWITCH] = int(use_oss)

class MonolithicVelocityPressureRansFormulation(RansFormulation):
    def __init__(self, model_part, settings):
        """Incompressible Variational-Multi-Scale Navier Stokes formulation

        This RansFormulation solves VELOCITY, and PRESSURE with Variational-Multi-Scale (VMS) formulated
        incompressible Navier-Stokes equation.

        This supports both steady and transient problems, transient with bossak time integration scheme.
        It uses wall functions at walls, therefore SLIP flag must be True on all walls (in both conditions and nodes)

        Args:
            model_part (Kratos.ModelPart): ModelPart to be used in the formulation.
            settings (Kratos.Parameters): Settings to be used in the formulation.
        """
        super().__init__(model_part, settings)

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
            "flow_solver_formulation": {
                "element_type": "vms",
                "use_orthogonal_subscales": false,
                "dynamic_tau": 0.01
            }
        }""")

        settings.ValidateAndAssignDefaults(default_settings)

        self.min_buffer_size = 2
        self.echo_level = settings["echo_level"].GetInt()

        self.flow_solver_formulation = StabilizedFormulation(settings["flow_solver_formulation"])
        self.flow_solver_formulation.SetProcessInfo(self.GetBaseModelPart())

        self.SetMaxCouplingIterations(1)

        Kratos.Logger.PrintInfo(self.__class__.__name__, "Construction of formulation finished.")

    def AddVariables(self):
        base_model_part = self.GetBaseModelPart()
        base_model_part.AddNodalSolutionStepVariable(Kratos.VELOCITY)
        base_model_part.AddNodalSolutionStepVariable(Kratos.ACCELERATION)
        base_model_part.AddNodalSolutionStepVariable(Kratos.MESH_VELOCITY)
        base_model_part.AddNodalSolutionStepVariable(Kratos.PRESSURE)
        base_model_part.AddNodalSolutionStepVariable(Kratos.IS_STRUCTURE)
        base_model_part.AddNodalSolutionStepVariable(Kratos.DISPLACEMENT)
        base_model_part.AddNodalSolutionStepVariable(Kratos.BODY_FORCE)
        base_model_part.AddNodalSolutionStepVariable(Kratos.NODAL_AREA)
        base_model_part.AddNodalSolutionStepVariable(Kratos.NODAL_H)
        base_model_part.AddNodalSolutionStepVariable(Kratos.ADVPROJ)
        base_model_part.AddNodalSolutionStepVariable(Kratos.DIVPROJ)
        base_model_part.AddNodalSolutionStepVariable(Kratos.REACTION)
        base_model_part.AddNodalSolutionStepVariable(Kratos.REACTION_WATER_PRESSURE)
        base_model_part.AddNodalSolutionStepVariable(Kratos.EXTERNAL_PRESSURE)
        base_model_part.AddNodalSolutionStepVariable(Kratos.NORMAL)
        base_model_part.AddNodalSolutionStepVariable(Kratos.Y_WALL)
        base_model_part.AddNodalSolutionStepVariable(KratosCFD.Q_VALUE)
        base_model_part.AddNodalSolutionStepVariable(KratosRANS.TURBULENT_KINETIC_ENERGY)

        if (self.flow_solver_formulation.element_name == "VMS"):
            base_model_part.AddNodalSolutionStepVariable(Kratos.DENSITY)
            base_model_part.AddNodalSolutionStepVariable(Kratos.VISCOSITY)

        Kratos.Logger.PrintInfo(self.__class__.__name__, "Added solution step variables.")

    def AddDofs(self):
        base_model_part = self.GetBaseModelPart()
        Kratos.VariableUtils().AddDof(Kratos.VELOCITY_X, Kratos.REACTION_X, base_model_part)
        Kratos.VariableUtils().AddDof(Kratos.VELOCITY_Y, Kratos.REACTION_Y, base_model_part)
        Kratos.VariableUtils().AddDof(Kratos.VELOCITY_Z, Kratos.REACTION_Z, base_model_part)
        Kratos.VariableUtils().AddDof(Kratos.PRESSURE, Kratos.REACTION_WATER_PRESSURE, base_model_part)

        Kratos.Logger.PrintInfo(self.__class__.__name__, "Added dofs.")

    def PrepareModelPart(self):
        self.monolithic_model_part = CreateRansFormulationModelPart(
            self.GetComputingModelPart(),
            self.__class__.__name__,
            self.GetDomainSize(),
            self.flow_solver_formulation.element_name,
            self.condition_name)
        Kratos.Logger.PrintInfo(self.__class__.__name__, "Created formulation model part.")

    def Initialize(self):
        model_part = self.GetBaseModelPart()
        CalculateNormalsOnConditions(model_part)

        process_info = model_part.ProcessInfo
        bossak_alpha = process_info[Kratos.BOSSAK_ALPHA]

        settings = self.GetParameters()

        if (self.IsPeriodic()):
            if (self.GetDomainSize() == 2):
                periodic_variables_list = [Kratos.VELOCITY_X, Kratos.VELOCITY_Y, Kratos.PRESSURE]
            else:
                periodic_variables_list = [Kratos.VELOCITY_X, Kratos.VELOCITY_Y, Kratos.VELOCITY_Z, Kratos.PRESSURE]
            InitializePeriodicConditions(model_part,
                                         self.monolithic_model_part,
                                         periodic_variables_list)

        conv_criteria_type = GetKratosObjectPrototype("MixedGenericCriteria")
        conv_criteria = conv_criteria_type(
            [(Kratos.VELOCITY, settings["relative_velocity_tolerance"].GetDouble(), settings["absolute_velocity_tolerance"].GetDouble()),
             (Kratos.PRESSURE, settings["relative_pressure_tolerance"].GetDouble(), settings["absolute_pressure_tolerance"].GetDouble())])

        if self.is_steady_simulation:
            scheme_type = GetKratosObjectPrototype("ResidualBasedSimpleSteadyScheme")
            scheme = scheme_type(
                settings["velocity_relaxation"].GetDouble(),
                settings["pressure_relaxation"].GetDouble(),
                self.GetDomainSize())
        else:
            scheme_type = GetKratosObjectPrototype("ResidualBasedPredictorCorrectorVelocityBossakSchemeTurbulent")
            scheme = scheme_type(
                bossak_alpha,
                settings["move_mesh_strategy"].GetInt(),
                self.GetDomainSize())

        linear_solver_factory = GetKratosObjectPrototype("LinearSolverFactory")
        linear_solver = linear_solver_factory(settings["linear_solver_settings"])

        builder_and_solver = CreateBlockBuilderAndSolver(
            linear_solver,
            self.IsPeriodic(),
            self.GetCommunicator())

        solver_type = GetKratosObjectPrototype("ResidualBasedNewtonRaphsonStrategy")
        self.solver = solver_type(
            self.monolithic_model_part,
            scheme,
            conv_criteria,
            builder_and_solver,
            settings["maximum_iterations"].GetInt(),
            settings["compute_reactions"].GetBool(),
            settings["reform_dofs_at_each_step"].GetBool(),
            settings["move_mesh_flag"].GetBool())

        self.solver.SetEchoLevel(self.echo_level)
        conv_criteria.SetEchoLevel(self.echo_level)

        super().Initialize()

        Kratos.Logger.PrintInfo(self.__class__.__name__, "Solver initialization finished.")

    def GetMinimumBufferSize(self):
        return self.min_buffer_size

    def SolveCouplingStep(self):
        if (self.IsBufferInitialized()):
            max_iterations = self.GetMaxCouplingIterations()
            for iteration in range(max_iterations):
                self.ExecuteBeforeCouplingSolveStep()
                self.solver.Predict()
                _ = self.solver.SolveSolutionStep()
                self.ExecuteAfterCouplingSolveStep()
                Kratos.Logger.PrintInfo(self.__class__.__name__, "Solved coupling iteration " + str(iteration + 1) + "/" + str(max_iterations) + ".")
                return True

        return False

    def InitializeSolutionStep(self):
        if (self.IsBufferInitialized()):
            super().InitializeSolutionStep()

    def FinalizeSolutionStep(self):
        if (self.IsBufferInitialized()):
            super().FinalizeSolutionStep()

    def SetTimeSchemeSettings(self, settings):
        if (settings.Has("scheme_type")):
            scheme_type = settings["scheme_type"].GetString()
            if (scheme_type == "steady"):
                self.is_steady_simulation = True
            elif (scheme_type == "bossak"):
                self.is_steady_simulation = False
                default_settings = Kratos.Parameters('''{
                    "scheme_type": "bossak",
                    "alpha_bossak": -0.3
                }''')
                settings.ValidateAndAssignDefaults(default_settings)
                self.GetBaseModelPart().ProcessInfo.SetValue(Kratos.BOSSAK_ALPHA, settings["alpha_bossak"].GetDouble())
            else:
                raise Exception("Only \"steady\" and \"bossak\" scheme types supported. [ scheme_type = \"" + scheme_type  + "\" ]")
        else:
            raise Exception("\"scheme_type\" is missing in time scheme settings")

    def SetConstants(self, settings):
        defaults = Kratos.Parameters('''{
            "von_karman": 0.41,
            "c_mu"      : 0.09
        }''')
        settings.ValidateAndAssignDefaults(defaults)

        # set constants
        von_karman = settings["von_karman"].GetDouble()

        process_info = self.GetBaseModelPart().ProcessInfo
        process_info.SetValue(KratosRANS.VON_KARMAN, von_karman)
        process_info.SetValue(KratosRANS.TURBULENCE_RANS_C_MU, settings["c_mu"].GetDouble())

    def SetWallFunctionSettings(self, settings):
        wall_function_region_type = "logarithmic_region_only"
        if (settings.Has("wall_function_region_type")):
            wall_function_region_type = settings["wall_function_region_type"].GetString()

        if (wall_function_region_type == "logarithmic_region_only"):
            self.condition_name = "RansVMSMonolithicKBasedWall"
        else:
            msg = "Unsupported wall function region type provided. [ wall_function_region_type = \"" + wall_function_region_type + "\" ]."
            msg += "Supported wall function region types are:\n"
            msg += "\tlogarithmic_region_only\n"
            raise Exception(msg)

    def GetStrategy(self):
        return self.solver

    def ElementHasNodalProperties(self):
        return self.flow_solver_formulation.element_has_nodal_properties

