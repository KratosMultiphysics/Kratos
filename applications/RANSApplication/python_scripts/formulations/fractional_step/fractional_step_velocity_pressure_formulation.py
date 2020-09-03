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
from KratosMultiphysics.RANSApplication import RansVariableUtilities
from KratosMultiphysics.RANSApplication import RansCalculationUtilities
from KratosMultiphysics.RANSApplication.formulations.utilities import CreateLinearSolver
from KratosMultiphysics.RANSApplication.formulations.utilities import CreateFormulationModelPart
from KratosMultiphysics.RANSApplication.formulations.utilities import CalculateNormalsOnConditions
from KratosMultiphysics.RANSApplication.formulations.utilities import IsBufferInitialized
from KratosMultiphysics.RANSApplication.formulations.utilities import GetConvergenceInfo

# case specific imports
if (IsDistributedRun() and CheckIfApplicationsAvailable("TrilinosApplication")):
    from KratosMultiphysics.FluidDynamicsApplication.TrilinosExtension import TrilinosFractionalStepSettingsPeriodic as periodic_solver_settings
    from KratosMultiphysics.FluidDynamicsApplication.TrilinosExtension import TrilinosFractionalStepSettings as solver_settings
    from KratosMultiphysics.RANSApplication.TrilinosExtension import MPIRansFractionalStepStrategy as solver_strategy
    from KratosMultiphysics.FluidDynamicsApplication.TrilinosExtension import TrilinosStrategyLabel as strategy_label
elif (not IsDistributedRun()):
    from KratosMultiphysics.FluidDynamicsApplication import FractionalStepSettingsPeriodic as periodic_solver_settings
    from KratosMultiphysics.FluidDynamicsApplication import FractionalStepSettings as solver_settings
    from KratosMultiphysics.RANSApplication import RansFractionalStepStrategy as solver_strategy
    from KratosMultiphysics.FluidDynamicsApplication import StrategyLabel as strategy_label
else:
    raise Exception("Distributed run requires TrilinosApplication")


def CreateFractionalStepSolverSettings(
                    is_periodic,
                    model_part,
                    domain_size,
                    time_order,
                    use_slip_conditions,
                    move_mesh_flag,
                    reform_dofs_at_each_step,
                    communicator):
    args = []
    if (IsDistributedRun()):
        args.append(communicator)
    args.extend([model_part, domain_size, time_order, use_slip_conditions, move_mesh_flag, reform_dofs_at_each_step])

    if (is_periodic):
        args.append(KratosCFD.PATCH_INDEX)
        return periodic_solver_settings(*args)
    else:
        return solver_settings(*args)

class FractionalStepVelocityPressureFormulation(Formulation):
    def __init__(self, model_part, settings):
        super(FractionalStepVelocityPressureFormulation, self).__init__(model_part, settings)

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
        self.element_has_nodal_properties = True
        self.fractional_step_model_part = None

        ## Construct the linear solvers
        self.pressure_linear_solver = CreateLinearSolver(self.settings["pressure_linear_solver_settings"])
        self.velocity_linear_solver = CreateLinearSolver(self.settings["velocity_linear_solver_settings"])

        self.velocity_tolerance = self.settings["velocity_tolerance"].GetDouble()
        self.pressure_tolerance = self.settings["pressure_tolerance"].GetDouble()
        self.echo_level = self.settings["echo_level"].GetInt()
        self.SetMaxCouplingIterations(1)

        self.compute_reactions = self.settings["compute_reactions"].GetBool()

        Kratos.Logger.PrintInfo(self.GetName(), "Construction of formulation finished.")

    def AddVariables(self):
        self.GetBaseModelPart().AddNodalSolutionStepVariable(Kratos.DENSITY)
        self.GetBaseModelPart().AddNodalSolutionStepVariable(Kratos.PRESSURE)
        self.GetBaseModelPart().AddNodalSolutionStepVariable(Kratos.VELOCITY)
        self.GetBaseModelPart().AddNodalSolutionStepVariable(Kratos.ACCELERATION)
        self.GetBaseModelPart().AddNodalSolutionStepVariable(Kratos.MESH_VELOCITY)
        self.GetBaseModelPart().AddNodalSolutionStepVariable(Kratos.BODY_FORCE)
        self.GetBaseModelPart().AddNodalSolutionStepVariable(Kratos.NODAL_H)
        self.GetBaseModelPart().AddNodalSolutionStepVariable(Kratos.REACTION)
        self.GetBaseModelPart().AddNodalSolutionStepVariable(Kratos.REACTION_WATER_PRESSURE)
        self.GetBaseModelPart().AddNodalSolutionStepVariable(Kratos.NORMAL)
        self.GetBaseModelPart().AddNodalSolutionStepVariable(Kratos.Y_WALL)
        self.GetBaseModelPart().AddNodalSolutionStepVariable(Kratos.EXTERNAL_PRESSURE)
        self.GetBaseModelPart().AddNodalSolutionStepVariable(Kratos.VISCOSITY)
        self.GetBaseModelPart().AddNodalSolutionStepVariable(Kratos.KINEMATIC_VISCOSITY)
        self.GetBaseModelPart().AddNodalSolutionStepVariable(KratosRANS.TURBULENT_KINETIC_ENERGY)
        self.GetBaseModelPart().AddNodalSolutionStepVariable(Kratos.FRACT_VEL)
        self.GetBaseModelPart().AddNodalSolutionStepVariable(Kratos.PRESSURE_OLD_IT)
        # The following are used for the calculation of projections
        self.GetBaseModelPart().AddNodalSolutionStepVariable(Kratos.NODAL_AREA)
        self.GetBaseModelPart().AddNodalSolutionStepVariable(Kratos.PRESS_PROJ)
        self.GetBaseModelPart().AddNodalSolutionStepVariable(Kratos.CONV_PROJ)
        self.GetBaseModelPart().AddNodalSolutionStepVariable(Kratos.DIVPROJ)
        self.GetBaseModelPart().AddNodalSolutionStepVariable(KratosCFD.Q_VALUE)

        Kratos.Logger.PrintInfo(self.GetName(), "Added solution step variables.")

    def AddDofs(self):
        Kratos.VariableUtils().AddDof(Kratos.VELOCITY_X, Kratos.REACTION_X,self.GetBaseModelPart())
        Kratos.VariableUtils().AddDof(Kratos.VELOCITY_Y, Kratos.REACTION_Y,self.GetBaseModelPart())
        Kratos.VariableUtils().AddDof(Kratos.VELOCITY_Z, Kratos.REACTION_Z,self.GetBaseModelPart())
        Kratos.VariableUtils().AddDof(Kratos.PRESSURE, Kratos.REACTION_WATER_PRESSURE,self.GetBaseModelPart())

        Kratos.Logger.PrintInfo(self.GetName(), "Added dofs.")

    def PrepareModelPart(self):
        if (self.is_steady_simulation):
            self.fractional_step_model_part = CreateFormulationModelPart(self,
                                                                        "RansFractionalStep",
                                                                        self.condition_name)
        else:
            self.fractional_step_model_part = CreateFormulationModelPart(self,
                                                                        "FractionalStep",
                                                                        self.condition_name)

        Kratos.Logger.PrintInfo(self.GetName(), "Created formulation model part.")

    def Initialize(self):
        model_part = self.GetBaseModelPart()
        CalculateNormalsOnConditions(model_part)

        process_info = model_part.ProcessInfo
        wall_model_part_name = process_info[KratosRANS.WALL_MODEL_PART_NAME]
        kappa = process_info[KratosRANS.WALL_VON_KARMAN]
        beta = process_info[KratosRANS.WALL_SMOOTHNESS_BETA]
        wall_function_update_process = KratosRANS.RansWallFunctionUpdateProcess(
                                                            model_part.GetModel(),
                                                            wall_model_part_name,
                                                            kappa,
                                                            beta,
                                                            self.echo_level)
        self.AddProcess(wall_function_update_process)

        if (self.IsPeriodic()):
            self.__InitializePeriodicConditions()

        self.solver_settings = CreateFractionalStepSolverSettings(
                                            self.IsPeriodic(),
                                            self.fractional_step_model_part,
                                            self.domain_size,
                                            self.settings["time_order"].GetInt(),
                                            True,
                                            self.GetMoveMeshFlag(),
                                            self.settings["reform_dofs_at_each_step"].GetBool(),
                                            self.GetCommunicator())

        self.solver_settings.SetEchoLevel(self.settings["echo_level"].GetInt())

        self.solver_settings.SetStrategy(strategy_label.Velocity,
                                         self.velocity_linear_solver,
                                         self.velocity_tolerance,
                                         self.settings["maximum_velocity_iterations"].GetInt())

        self.solver_settings.SetStrategy(strategy_label.Pressure,
                                         self.pressure_linear_solver,
                                         self.pressure_tolerance,
                                         self.settings["maximum_pressure_iterations"].GetInt())


        if self.IsPeriodic():
            self.solver = solver_strategy(self.fractional_step_model_part,
                                          self.solver_settings,
                                          self.settings["predictor_corrector"].GetBool(),
                                          False,
                                          KratosCFD.PATCH_INDEX)
        else:
            self.solver = solver_strategy(self.fractional_step_model_part,
                                          self.solver_settings,
                                          self.settings["predictor_corrector"].GetBool(),
                                          False)

        self.GetBaseModelPart().ProcessInfo.SetValue(Kratos.DYNAMIC_TAU, self.settings["dynamic_tau"].GetDouble())
        self.GetBaseModelPart().ProcessInfo.SetValue(Kratos.OSS_SWITCH, self.settings["oss_switch"].GetInt())

        super(FractionalStepVelocityPressureFormulation, self).Initialize()
        self.solver.Initialize()

        Kratos.Logger.PrintInfo(self.GetName(), "Solver initialization finished.")

    def GetMinimumBufferSize(self):
        return self.min_buffer_size

    def Finalize(self):
        self.solver.Clear()
        super(FractionalStepVelocityPressureFormulation, self).Finalize()

    def IsConverged(self):
        if (hasattr(self, "is_converged")):
            if (self.is_steady_simulation):
                settings = self.settings["steady_convergence_settings"]
                formulation_converged = self.__CheckTransientConvergence(Kratos.VELOCITY, settings["velocity_tolerances"])
                self.is_converged = self.is_converged and formulation_converged

                formulation_converged = self.__CheckTransientConvergence(Kratos.PRESSURE, settings["pressure_tolerances"])
                self.is_converged = self.is_converged and formulation_converged

            return self.is_converged
        return False

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
            super(FractionalStepVelocityPressureFormulation, self).InitializeSolutionStep()
            self.solver.InitializeSolutionStep()

    def FinializeSolutionStep(self):
        if (IsBufferInitialized(self)):
            self.solver.FinializeSolutionStep()
            super(FractionalStepVelocityPressureFormulation, self).FinializeSolutionStep()

    def Check(self):
        super(FractionalStepVelocityPressureFormulation, self).Check()
        self.solver.Check()

    def Clear(self):
        self.solver.Clear()
        super(FractionalStepVelocityPressureFormulation, self).Clear()

    def SetTimeSchemeSettings(self, settings):
        if (settings.Has("scheme_type")):
            scheme_type = settings["scheme_type"].GetString()
            if (scheme_type == "steady"):
                self.is_steady_simulation = True
                default_settings = Kratos.Parameters('''{
                    "scheme_type": "steady",
                    "pressure_coefficient": 0.5
                }''')
                settings.ValidateAndAssignDefaults(default_settings)
                self.GetBaseModelPart().ProcessInfo.SetValue(Kratos.PRESSURE_COEFFICIENT, settings["pressure_coefficient"].GetDouble())
            elif (scheme_type == "transient"):
                self.is_steady_simulation = False
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

    def GetStrategy(self):
        return self.solver

    def GetModelPart(self):
        return self.fractional_step_model_part

    def SetWallFunctionSettings(self, settings):
        wall_function_region_type = "logarithmic_region_only"
        if (settings.Has("wall_function_region_type")):
            wall_function_region_type = settings[
                "wall_function_region_type"].GetString()

        if (wall_function_region_type == "logarithmic_region_only"):
            self.condition_name = "RansFractionalStepKBasedWall"
        else:
            msg = "Unsupported wall function region type provided. [ wall_function_region_type = \"" + wall_function_region_type + "\" ]."
            msg += "Supported wall function region types are:\n"
            msg += "\tlogarithmic_region_only\n"
            raise Exception(msg)

    def __CheckTransientConvergence(self, variable, settings):
        relative_error, absolute_error = RansVariableUtilities.CalculateTransientVariableConvergence(self.GetBaseModelPart(), variable)
        relative_tolerance = settings["relative_tolerance"].GetDouble()
        absolute_tolerance = settings["absolute_tolerance"].GetDouble()

        info = GetConvergenceInfo(variable, relative_error, relative_tolerance, absolute_error, absolute_tolerance)
        print("+++++++")
        print(variable)
        print(relative_error)
        print(absolute_error)
        Kratos.Logger.PrintInfo(self.GetName(), info)

        return (relative_error <= relative_tolerance or absolute_error <= absolute_tolerance)

    def __InitializePeriodicConditions(self):
        base_model_part = self.GetBaseModelPart()
        model_part = self.fractional_step_model_part
        properties = model_part.CreateNewProperties(
            model_part.NumberOfProperties() + 1)

        periodic_condition_name = "FSPeriodicCondition{0:d}D".format(self.domain_size)

        current_number_of_conditions = model_part.NumberOfConditions()
        periodic_conditions_count = 0
        for condition in base_model_part.Conditions:
            if condition.Is(Kratos.PERIODIC):
                periodic_conditions_count += 1
                node_id_list = [node.Id for node in condition.GetNodes()]
                periodic_condition = model_part.CreateNewCondition(
                    periodic_condition_name, current_number_of_conditions + periodic_conditions_count, node_id_list, properties)
                periodic_condition.Set(Kratos.PERIODIC)

        Kratos.Logger.PrintInfo(self.GetName(), "Created {0:d} periodic conditions using {1:s}".format(periodic_conditions_count, periodic_condition_name))
