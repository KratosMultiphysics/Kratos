from pathlib import Path

# Importing the Kratos Library
import KratosMultiphysics as Kratos

# Import applications
import KratosMultiphysics.FluidDynamicsApplication as KratosCFD
import KratosMultiphysics.RANSApplication as KratosRANS
from KratosMultiphysics.kratos_utilities import IssueDeprecationWarning

# import formulation interface
from KratosMultiphysics.RANSApplication.formulations.rans_formulation import RansFormulation

# import utilities
from KratosMultiphysics.RANSApplication.formulations.utilities import SolveProblem
from KratosMultiphysics.RANSApplication.formulations.utilities import CreateBlockBuilderAndSolver
from KratosMultiphysics.RANSApplication.formulations.utilities import CreateRansFormulationModelPart
from KratosMultiphysics.RANSApplication.formulations.utilities import CalculateNormalsOnConditions
from KratosMultiphysics.RANSApplication.formulations.utilities import InitializePeriodicConditions
from KratosMultiphysics.RANSApplication.formulations.utilities import GetKratosObjectPrototype
from KratosMultiphysics.RANSApplication.formulations.utilities import ExecutionScope
from KratosMultiphysics.RANSApplication.formulations.utilities import AddWallPropertiesUpdateProcess

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
    def __init__(self, model_part, settings, deprecated_settings_dict):
        """Incompressible Variational-Multi-Scale Navier Stokes formulation

        This RansFormulation solves VELOCITY, and PRESSURE with Variational-Multi-Scale (VMS) formulated
        incompressible Navier-Stokes equation.

        This supports both steady and transient problems, transient with bossak time integration scheme.
        It uses wall functions at walls, therefore SLIP flag must be True on all walls (in both conditions and nodes)

        Args:
            model_part (Kratos.ModelPart): ModelPart to be used in the formulation.
            settings (Kratos.Parameters): Settings to be used in the formulation.
        """
        self.BackwardCompatibilityHelper(settings, deprecated_settings_dict)
        super().__init__(model_part, settings)

        settings.ValidateAndAssignDefaults(self.GetDefaultParameters())

        self.min_buffer_size = 2
        self.echo_level = settings["echo_level"].GetInt()

        self.flow_solver_formulation = StabilizedFormulation(settings["flow_solver_formulation"])
        self.flow_solver_formulation.SetProcessInfo(self.GetBaseModelPart())

        self.SetMaxCouplingIterations(1)

        Kratos.Logger.PrintInfo(self.__class__.__name__, "Construction of formulation finished.")

    def GetDefaultParameters(self):
        return Kratos.Parameters("""
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
            },
            "use_frozen_turbulence" : false,
            "additional_frozen_turbulence_constants": {
                "a1"  : 0.31
            },
            "wall_function_settings": {
                "wall_friction_velocity_calculation_method": "turbulent_kinetic_energy_based",
                "wall_function_region_type": "logarithmic_region_only"
            }
        }""")

    def BackwardCompatibilityHelper(self, settings, deprecated_settings_dict):
        if "wall_function_settings" in deprecated_settings_dict.keys():
            if settings.Has("wall_function_settings"):
                Kratos.Logger.PrintWarning(self.__class__.__name__, "Found \"wall_function_settings\" in deprecated settings as well as in formulation settings. Continuing with formulation based settings.")
            else:
                settings.AddValue("wall_function_settings", deprecated_settings_dict["wall_function_settings"].Clone())

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

        if (self.GetParameters()["use_frozen_turbulence"].GetBool()):
            Kratos.Logger.PrintInfo(self.__class__.__name__, "Frozen turbulence assumption is used. Adding respective variables.")
            base_model_part.AddNodalSolutionStepVariable(KratosRANS.TURBULENT_ENERGY_DISSIPATION_RATE)
            base_model_part.AddNodalSolutionStepVariable(KratosRANS.TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE)
            base_model_part.AddNodalSolutionStepVariable(Kratos.DISTANCE)

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

        for constraint in model_part.MasterSlaveConstraints:
            if (constraint.GetSlaveDofsVector()[0].GetVariable() == Kratos.VELOCITY_X or
                constraint.GetSlaveDofsVector()[0].GetVariable() == Kratos.VELOCITY_Y or
                constraint.GetSlaveDofsVector()[0].GetVariable() == Kratos.VELOCITY_Z or
                constraint.GetSlaveDofsVector()[0].GetVariable() == Kratos.PRESSURE):
                self.monolithic_model_part.AddMasterSlaveConstraint(constraint)

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
            self.scheme = scheme_type(
                settings["velocity_relaxation"].GetDouble(),
                settings["pressure_relaxation"].GetDouble(),
                self.GetDomainSize())
        else:
            scheme_type = GetKratosObjectPrototype("ResidualBasedPredictorCorrectorVelocityBossakSchemeTurbulent")
            self.scheme = scheme_type(
                bossak_alpha,
                settings["move_mesh_strategy"].GetInt(),
                self.GetDomainSize())

        linear_solver_factory = GetKratosObjectPrototype("LinearSolverFactory")
        linear_solver = linear_solver_factory(settings["linear_solver_settings"])

        builder_and_solver = CreateBlockBuilderAndSolver(
            linear_solver,
            self.IsPeriodic(),
            self.GetCommunicator())

        self.compute_reactions_using_scheme = False
        if settings["compute_reactions"].GetBool():
            if hasattr(self.scheme, "CalculateReactions"):
                self.compute_reactions_using_scheme = True
            else:
                Kratos.Logger.PrintWarning(
                    self.__class__.__name__, "Using reaction computation from builder and solver because {:s} has not implemented CalculateReactions method. This will give wrong REACTION on nodes if wall functions are used.".format(self.scheme.__class__.__name__))

        solver_type = GetKratosObjectPrototype("ResidualBasedNewtonRaphsonStrategy")
        self.solver = solver_type(
            self.monolithic_model_part,
            self.scheme,
            conv_criteria,
            builder_and_solver,
            settings["maximum_iterations"].GetInt(),
            settings["compute_reactions"].GetBool() and not self.compute_reactions_using_scheme,
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
                _ = self.solver.SolveSolutionStep()
                if self.compute_reactions_using_scheme:
                    self.scheme.CalculateReactions(self.GetModelPart())
                    Kratos.Logger.PrintInfo(self.__class__.__name__, "{:s} is used to calculate reactions in {:s}.".format(self.scheme.__class__.__name__, self.GetModelPart().FullName()))
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

        if (self.GetParameters()["use_frozen_turbulence"].GetBool()):
            process_info.SetValue(KratosRANS.TURBULENCE_RANS_A1, self.GetParameters()["additional_frozen_turbulence_constants"]["a1"].GetDouble())
            Kratos.Logger.PrintInfo(self.__class__.__name__, "Frozen turbulence assumption is used. Added respective constants.")

    def SetWallFunctionSettings(self):
        wall_function_settings = self.GetParameters()["wall_function_settings"]
        wall_function_settings.ValidateAndAssignDefaults(self.GetDefaultParameters()["wall_function_settings"])
        self.condition_name = self.GetConditionNamePrefix()

        if (self.condition_name != ""):
            wall_function_region_type = wall_function_settings["wall_function_region_type"].GetString()
            wall_friction_velocity_calculation_method = wall_function_settings["wall_friction_velocity_calculation_method"].GetString()

            if (wall_function_region_type == "logarithmic_region_only"):
                if (wall_friction_velocity_calculation_method == "velocity_based"):
                    self.condition_name = self.condition_name + "UBasedWall"
                elif (wall_friction_velocity_calculation_method ==
                    "turbulent_kinetic_energy_based"):
                    self.condition_name = self.condition_name + "KBasedWall"
                else:
                    msg = "Unsupported wall friction velocity calculation method. [ wall_friction_velocity_calculation_method = \"" + wall_friction_velocity_calculation_method + "\" ].\n"
                    msg += "Supported methods are:\n"
                    msg += "\tvelocity_based\n"
                    msg += "\tturbulent_kinetic_energy_based\n"
                    raise Exception(msg)
            else:
                msg = "Unsupported wall function region type provided. [ wall_function_region_type = \"" + wall_function_region_type + "\" ]."
                msg += "Supported wall function region types are:\n"
                msg += "\tlogarithmic_region_only\n"
                raise Exception(msg)

        AddWallPropertiesUpdateProcess(self, wall_function_settings)

    def GetStrategy(self):
        if (hasattr(self, "solver")):
            return self.solver
        else:
            return None

    def ElementHasNodalProperties(self):
        return self.flow_solver_formulation.element_has_nodal_properties

    def GetElementNames(self):
        return [self.flow_solver_formulation.element_name]

    def GetConditionNames(self):
        return [self.condition_name]

    def GetSolvingVariables(self):
        if (self.GetDomainSize() == 2):
            return [Kratos.VELOCITY_X, Kratos.VELOCITY_Y, Kratos.PRESSURE]
        else:
            return [Kratos.VELOCITY_X, Kratos.VELOCITY_Y, Kratos.VELOCITY_Z, Kratos.PRESSURE]


    def ComputeTransientResponseFunctionInterpolationError(self, settings, amr_output_path, model_import_settings):
        # in here we need to calculate response function interpolation error
        default_settings = Kratos.Parameters("""{
            "primal_transient_project_parameters_file_name" : "PLEASE_SPECIFY_TRANSIENT_PRIMAL_PROJECT_PARAMETERS_FILE",
            "adjoint_transient_project_parameters_file_name": "PLEASE_SPECIFY_TRANSIENT_ADJOINT_PROJECT_PARAMETERS_FILE",
            "number_of_transient_steps_to_consider"         : 50
        }""")

        settings.ValidateAndAssignDefaults(default_settings)

        current_model_part = self.GetBaseModelPart()
        current_model = current_model_part.GetModel()
        process_info = current_model_part.ProcessInfo
        current_time = process_info[Kratos.TIME]
        current_step = process_info[Kratos.STEP]
        current_dt   = process_info[Kratos.DELTA_TIME]

        execution_path = Path(amr_output_path) / "step_{:d}".format(current_step)
        execution_path.mkdir(exist_ok=True, parents=True)

        # run the forward primal solution for user defined steps in the transient problem
        # while calculating time averaged quantities
        time_steps = settings["number_of_transient_steps_to_consider"].GetInt()
        start_time = current_time
        end_time = current_time + current_dt * time_steps

        Kratos.Logger.PrintInfo(self.__class__.__name__, "Solving forward primal solution for transient response function interpolation error calculation...")

        # open primal parameters json file
        primal_transient_project_parameters_file_name = settings["primal_transient_project_parameters_file_name"].GetString()
        with open(primal_transient_project_parameters_file_name, "r") as file_input:
            primal_parameters = Kratos.Parameters(file_input.read().replace("<error_computation_step>", str(current_step)))

        # set start time and end time
        primal_parameters["problem_data"]["start_time"].SetDouble(start_time)
        primal_parameters["problem_data"]["end_time"].SetDouble(end_time)
        primal_parameters["solver_settings"]["time_stepping"]["time_step"].SetDouble(current_dt)

        # solve primal problem
        with ExecutionScope(execution_path):
            from KratosMultiphysics.RANSApplication.rans_analysis import RANSAnalysis
            # write existing model part
            Kratos.ModelPartIO(primal_parameters["solver_settings"]["model_import_settings"]["input_filename"].GetString(), Kratos.IO.WRITE | Kratos.IO.MESH_ONLY).WriteModelPart(self.GetBaseModelPart())
            primal_model, primal_simulation = SolveProblem(RANSAnalysis, primal_parameters, primal_transient_project_parameters_file_name[:-5])

            # copy time averaged quantities from the primal_simulation
            time_averaged_variable_data = [
                ("TIME_AVERAGED_{:s}".format(var.Name()), False, 0, "TIME_AVERAGED_{:s}".format(var.Name()), False, 0) for var in self.GetSolvingVariables()
            ]
            KratosRANS.RansVariableDataTransferProcess(
                primal_model,
                current_model,
                primal_simulation._GetSolver().GetComputingModelPart().FullName(),
                current_model_part.FullName(),
                ["execute"],
                time_averaged_variable_data,
                self.echo_level).Execute()

            # free the memory consumed by primal analysis
            del primal_simulation
            del primal_model

        # now run the backward adjoint steady problem
        Kratos.Logger.PrintInfo(self.__class__.__name__, "Solving adjoint solution for transient response function interpolation error calculation...")

        # open adjoint parameters json file
        adjoint_transient_project_parameters_file_name = settings["adjoint_transient_project_parameters_file_name"].GetString()
        with open(adjoint_transient_project_parameters_file_name, "r") as file_input:
            adjoint_parameters = Kratos.Parameters(file_input.read().replace("<error_computation_step>", str(current_step)))

        # set start time and end time
        adjoint_parameters["problem_data"]["start_time"].SetDouble(start_time + current_dt)
        adjoint_parameters["problem_data"]["end_time"].SetDouble(end_time)

        # solve adjoint problem
        with ExecutionScope(execution_path):
            from KratosMultiphysics.RANSApplication.adjoint_rans_analysis import AdjointRANSAnalysis
            adjoint_model, adjoint_simulation = SolveProblem(AdjointRANSAnalysis, adjoint_parameters, adjoint_transient_project_parameters_file_name[:-5])

            # now transfer RESPONSE_FUNCTION_INTERPOLATION_ERROR data to the current model part
            KratosRANS.RansVariableDataTransferProcess(
                adjoint_model,
                current_model,
                adjoint_simulation._GetSolver().GetComputingModelPart().FullName(),
                current_model_part.FullName(),
                ["execute"],
                [("RESPONSE_FUNCTION_INTERPOLATION_ERROR", False, 0, "RESPONSE_FUNCTION_INTERPOLATION_ERROR", False, 0)],
                self.echo_level).Execute()

            # free the memory consumed by primal analysis
            del adjoint_simulation
            del adjoint_model

        Kratos.Logger.PrintInfo(self.__class__.__name__, "Computed response based interpolation errors for adaptive mesh refinement.")
    def GetConditionNamePrefix(self):
        return "RansVMSMonolithic"