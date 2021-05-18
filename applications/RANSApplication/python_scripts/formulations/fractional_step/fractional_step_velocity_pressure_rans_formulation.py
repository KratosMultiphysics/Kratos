# Importing the Kratos Library
import KratosMultiphysics as Kratos
from KratosMultiphysics import IsDistributedRun

# Import applications
import KratosMultiphysics.FluidDynamicsApplication as KratosCFD
import KratosMultiphysics.RANSApplication as KratosRANS

# import formulation interface
from KratosMultiphysics.RANSApplication.formulations.rans_formulation import RansFormulation

# import utilities
from KratosMultiphysics.RANSApplication import RansVariableUtilities
from KratosMultiphysics.RANSApplication.formulations.utilities import CreateRansFormulationModelPart
from KratosMultiphysics.RANSApplication.formulations.utilities import CalculateNormalsOnConditions
from KratosMultiphysics.RANSApplication.formulations.utilities import GetConvergenceInfo
from KratosMultiphysics.RANSApplication.formulations.utilities import GetKratosObjectPrototype
from KratosMultiphysics.RANSApplication.formulations.utilities import InitializePeriodicConditions


class FractionalStepVelocityPressureRansFormulation(RansFormulation):
    def __init__(self, model_part, settings):
        """Incompressible Fractional Step Navier Stokes formulation

        This RansFormulation solves VELOCITY, and PRESSURE with Fractional Step formulated
        incompressible Navier-Stokes equation.

        This supports both steady and transient problems, transient with time integration scheme.
        It uses wall functions at walls, therefore SLIP flag must be True on all walls (in both conditions and nodes)

        Args:
            model_part (Kratos.ModelPart): ModelPart to be used in the formulation.
            settings (Kratos.Parameters): Settings to be used in the formulation.
        """
        super().__init__(model_part, settings)

        ##settings string in json format
        default_settings = Kratos.Parameters("""
        {
            "formulation_name": "fractional_step",
            "predictor_corrector": false,
            "maximum_velocity_iterations": 3,
            "maximum_pressure_iterations": 3,
            "velocity_tolerance": 1e-3,
            "pressure_tolerance": 1e-2,
            "dynamic_tau": 0.01,
            "oss_switch": 0,
            "echo_level": 0,
            "time_order": 2,
            "compute_reactions": false,
            "reform_dofs_at_each_step": false,
            "pressure_linear_solver_settings":  {
                "solver_type"                    : "amgcl",
                "max_iteration"                  : 200,
                "tolerance"                      : 1e-6,
                "provide_coordinates"            : false,
                "smoother_type"                  : "ilu0",
                "krylov_type"                    : "cg",
                "gmres_krylov_space_dimension"   : 100,
                "use_block_matrices_if_possible" : false,
                "coarsening_type"                : "aggregation",
                "scaling"                        : true,
                "verbosity"                      : 0
            },
            "velocity_linear_solver_settings": {
                "solver_type"                    : "amgcl",
                "max_iteration"                  : 200,
                "tolerance"                      : 1e-6,
                "provide_coordinates"            : false,
                "smoother_type"                  : "ilu0",
                "krylov_type"                    : "lgmres",
                "gmres_krylov_space_dimension"   : 100,
                "use_block_matrices_if_possible" : false,
                "coarsening_type"                : "aggregation",
                "scaling"                        : true,
                "verbosity"                      : 0
            },
            "steady_convergence_settings":
            {
                "velocity_tolerances": {
                    "relative_tolerance": 1e-3,
                    "absolute_tolerance": 1e-5
                },
                "pressure_tolerances": {
                    "relative_tolerance": 1e-3,
                    "absolute_tolerance": 1e-5
                }
            }
        }""")

        settings.ValidateAndAssignDefaults(default_settings)

        self.min_buffer_size = 3

        self.echo_level = settings["echo_level"].GetInt()
        self.SetMaxCouplingIterations(1)

        Kratos.Logger.PrintInfo(self.__class__.__name__, "Construction of formulation finished.")

    def AddVariables(self):
        base_model_part = self.GetBaseModelPart()
        base_model_part.AddNodalSolutionStepVariable(Kratos.DENSITY)
        base_model_part.AddNodalSolutionStepVariable(Kratos.PRESSURE)
        base_model_part.AddNodalSolutionStepVariable(Kratos.VELOCITY)
        base_model_part.AddNodalSolutionStepVariable(Kratos.ACCELERATION)
        base_model_part.AddNodalSolutionStepVariable(Kratos.MESH_VELOCITY)
        base_model_part.AddNodalSolutionStepVariable(Kratos.BODY_FORCE)
        base_model_part.AddNodalSolutionStepVariable(Kratos.NODAL_H)
        base_model_part.AddNodalSolutionStepVariable(Kratos.REACTION)
        base_model_part.AddNodalSolutionStepVariable(Kratos.REACTION_WATER_PRESSURE)
        base_model_part.AddNodalSolutionStepVariable(Kratos.NORMAL)
        base_model_part.AddNodalSolutionStepVariable(Kratos.Y_WALL)
        base_model_part.AddNodalSolutionStepVariable(Kratos.EXTERNAL_PRESSURE)
        base_model_part.AddNodalSolutionStepVariable(Kratos.VISCOSITY)
        base_model_part.AddNodalSolutionStepVariable(KratosRANS.TURBULENT_KINETIC_ENERGY)
        base_model_part.AddNodalSolutionStepVariable(Kratos.FRACT_VEL)
        base_model_part.AddNodalSolutionStepVariable(Kratos.PRESSURE_OLD_IT)
        # The following are used for the calculation of projections
        base_model_part.AddNodalSolutionStepVariable(Kratos.NODAL_AREA)
        base_model_part.AddNodalSolutionStepVariable(Kratos.PRESS_PROJ)
        base_model_part.AddNodalSolutionStepVariable(Kratos.CONV_PROJ)
        base_model_part.AddNodalSolutionStepVariable(Kratos.DIVPROJ)
        base_model_part.AddNodalSolutionStepVariable(KratosCFD.Q_VALUE)

        Kratos.Logger.PrintInfo(self.__class__.__name__, "Added solution step variables.")

    def AddDofs(self):
        base_model_part = self.GetBaseModelPart()
        Kratos.VariableUtils().AddDof(Kratos.VELOCITY_X, Kratos.REACTION_X, base_model_part)
        Kratos.VariableUtils().AddDof(Kratos.VELOCITY_Y, Kratos.REACTION_Y, base_model_part)
        Kratos.VariableUtils().AddDof(Kratos.VELOCITY_Z, Kratos.REACTION_Z, base_model_part)
        Kratos.VariableUtils().AddDof(Kratos.PRESSURE, Kratos.REACTION_WATER_PRESSURE, base_model_part)

        Kratos.Logger.PrintInfo(self.__class__.__name__, "Added dofs.")

    def PrepareModelPart(self):
        self.fractional_step_model_part = CreateRansFormulationModelPart(
            self.GetComputingModelPart(),
            self.__class__.__name__,
            self.GetDomainSize(),
            "FractionalStep",
            self.condition_name)

        Kratos.Logger.PrintInfo(self.__class__.__name__, "Created formulation model part.")

    def Initialize(self):
        model_part = self.GetBaseModelPart()
        CalculateNormalsOnConditions(model_part)

        process_info = model_part.ProcessInfo
        wall_model_part_name = process_info[KratosRANS.WALL_MODEL_PART_NAME]
        wall_function_update_process = KratosRANS.RansWallFunctionUpdateProcess(
            model_part.GetModel(),
            wall_model_part_name,
            self.echo_level)

        self.AddProcess(wall_function_update_process)

        if (self.IsPeriodic()):
            InitializePeriodicConditions(
                model_part,
                self.fractional_step_model_part,
                None,
                "FSPeriodicCondition{0:d}D".format(self.GetDomainSize()))

        settings = self.GetParameters()

        self.solver_settings = self._CreateSolverSettings(
            self.fractional_step_model_part,
            self.GetDomainSize(),
            settings["time_order"].GetInt(),
            True,
            self.GetMoveMeshFlag(),
            settings["reform_dofs_at_each_step"].GetBool())

        self.solver_settings.SetEchoLevel(self.echo_level)

        ## Construct the linear solvers
        linear_solver_factory = GetKratosObjectPrototype("LinearSolverFactory")
        pressure_linear_solver = linear_solver_factory(settings["pressure_linear_solver_settings"])
        velocity_linear_solver = linear_solver_factory(settings["velocity_linear_solver_settings"])

        strategy_label_type = GetKratosObjectPrototype("StrategyLabel")

        self.solver_settings.SetStrategy(
            strategy_label_type.Velocity,
            velocity_linear_solver,
            settings["velocity_tolerance"].GetDouble(),
            settings["maximum_velocity_iterations"].GetInt())
        self.solver_settings.SetStrategy(
            strategy_label_type.Pressure,
            pressure_linear_solver,
            settings["pressure_tolerance"].GetDouble(),
            settings["maximum_pressure_iterations"].GetInt())

        solver_type = GetKratosObjectPrototype("FractionalStepStrategy")
        if self.IsPeriodic():
            self.solver = solver_type(
                self.fractional_step_model_part,
                self.solver_settings,
                settings["predictor_corrector"].GetBool(),
                False,
                KratosCFD.PATCH_INDEX)
        else:
            self.solver = solver_type(
                self.fractional_step_model_part,
                self.solver_settings,
                settings["predictor_corrector"].GetBool(),
                False)

        process_info.SetValue(Kratos.DYNAMIC_TAU, settings["dynamic_tau"].GetDouble())
        process_info.SetValue(Kratos.OSS_SWITCH, settings["oss_switch"].GetInt())

        super().Initialize()

        Kratos.Logger.PrintInfo(self.__class__.__name__, "Solver initialization finished.")

    def GetMinimumBufferSize(self):
        return self.min_buffer_size

    def IsConverged(self):
        if (hasattr(self, "is_converged")):
            if (self.is_steady_simulation):
                settings = self.GetParameters()["steady_convergence_settings"]
                formulation_converged = self._CheckTransientConvergence(
                    Kratos.VELOCITY,
                    settings["velocity_tolerances"])
                self.is_converged = self.is_converged and formulation_converged

                formulation_converged = self._CheckTransientConvergence(
                    Kratos.PRESSURE,
                    settings["pressure_tolerances"])
                self.is_converged = self.is_converged and formulation_converged

            return self.is_converged
        return False

    def SolveCouplingStep(self):
        if (self.IsBufferInitialized()):
            max_iterations = self.GetMaxCouplingIterations()
            for iteration in range(max_iterations):
                self.ExecuteBeforeCouplingSolveStep()
                self.solver.Predict()
                self.is_converged = self.solver.SolveSolutionStep()
                self.ExecuteAfterCouplingSolveStep()
                Kratos.Logger.PrintInfo(self.__class__.__name__, "Solved coupling iteration " + str(iteration + 1) + "/" + str(max_iterations) + ".")
                return True

        return False

    def SetTimeSchemeSettings(self, settings):
        if (settings.Has("scheme_type")):
            scheme_type = settings["scheme_type"].GetString()
            if (scheme_type == "steady"):
                self.is_steady_simulation = True
                default_settings = Kratos.Parameters('''{
                    "scheme_type": "steady",
                    "pressure_gradient_relaxation_factor": 0.5
                }''')
                settings.ValidateAndAssignDefaults(default_settings)
                self.GetBaseModelPart().ProcessInfo.SetValue(KratosCFD.FS_PRESSURE_GRADIENT_RELAXATION_FACTOR, settings["pressure_gradient_relaxation_factor"].GetDouble())
            elif (scheme_type == "bdf2"):
                self.is_steady_simulation = False
            else:
                raise Exception("Only \"steady\" and \"bdf2\" scheme types supported. [ scheme_type = \"" + scheme_type  + "\" ]")
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

    def GetStrategy(self):
        return self.solver

    def SetWallFunctionSettings(self, settings):
        wall_function_region_type = "logarithmic_region_only"
        if (settings.Has("wall_function_region_type")):
            wall_function_region_type = settings["wall_function_region_type"].GetString()

        if (wall_function_region_type == "logarithmic_region_only"):
            self.condition_name = "RansFractionalStepKBasedWall"
        else:
            msg = "Unsupported wall function region type provided. [ wall_function_region_type = \"" + wall_function_region_type + "\" ]."
            msg += "Supported wall function region types are:\n"
            msg += "\tlogarithmic_region_only\n"
            raise Exception(msg)

    def _CheckTransientConvergence(self, variable, settings):
        relative_error, absolute_error = RansVariableUtilities.CalculateTransientVariableConvergence(
            self.GetBaseModelPart(),
            variable)
        relative_tolerance = settings["relative_tolerance"].GetDouble()
        absolute_tolerance = settings["absolute_tolerance"].GetDouble()

        info = GetConvergenceInfo(
            variable,
            relative_error,
            relative_tolerance,
            absolute_error,
            absolute_tolerance)
        Kratos.Logger.PrintInfo(self.__class__.__name__, info)

        return (relative_error <= relative_tolerance or absolute_error <= absolute_tolerance)

    def _CreateSolverSettings(self, *args):
        if (self.IsPeriodic()):
            solver_settings_type = GetKratosObjectPrototype("FractionalStepSettingsPeriodic")
        else:
            solver_settings_type = GetKratosObjectPrototype("FractionalStepSettings")

        if (IsDistributedRun()):
            return solver_settings_type(self.GetCommunicator(), *args)
        else:
            return solver_settings_type(*args)
