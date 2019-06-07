import KratosMultiphysics as Kratos
import KratosMultiphysics.RANSModellingApplication as KratosRANS
import math

from KratosMultiphysics.kratos_utilities import CheckIfApplicationsAvailable

if CheckIfApplicationsAvailable("FluidDynamicsApplication"):
    import KratosMultiphysics.FluidDynamicsApplication as KratosCFD
    from turbulence_model_configuration import TurbulenceModelConfiguration
else:
    msg = "RANSModellingApplication requires FluidDynamicsApplication which is not found."
    msg += " Please install/compile it and try again."
    raise Exception(msg)


class TurbulenceEddyViscosityModelConfiguration(TurbulenceModelConfiguration):
    def __init__(self, model, settings):
        self._validate_settings_in_baseclass = True  # To be removed eventually

        super(TurbulenceEddyViscosityModelConfiguration, self).__init__(
            model, settings)

        # self.mesh_moving = self.settings["mesh_moving"].GetBool()
        self.distance_calculation_process = None
        self.y_plus_model_process = None
        self.wall_velocity_model_process = None
        self.turbulence_model_process = None
        self.strategies_list = []
        self.nu_t_min = self.settings["turbulent_viscosity_min"].GetDouble()
        self.nu_t_max = self.settings["turbulent_viscosity_max"].GetDouble()

        # TODO: Implement stuff for mesh_moving

        self.is_computing_solution = False

        self.model_elements_list = []
        self.model_conditions_list = []
        self.model_parts_list = []


    def GetDefaultSettings(self):
        return Kratos.Parameters(r'''{
            "model_type"            : "",
            "velocity_pressure_relaxation_factor": 0.2,
            "model_settings"        : {},
            "wall_distance_calculation"  : {
                "model_type"     : "",
                "model_settings" : {}
            },
            "y_plus_model"      : {
                "model_type"     : "",
                "model_settings" : {}
            },
            "wall_velocity_calculation":{
                "model_type"     : "",
                "model_settings" : {}
            },
            "mesh_moving"             : false,
            "echo_level"              : 0,
            "turbulent_viscosity_min" : 1e-12,
            "turbulent_viscosity_max" : 1e+2
        }''')

    def AddVariables(self):
        # adding variables required by rans eddy viscosity models
        self.fluid_model_part.AddNodalSolutionStepVariable(Kratos.DISTANCE)
        self.fluid_model_part.AddNodalSolutionStepVariable(
            Kratos.FLAG_VARIABLE)
        self.fluid_model_part.AddNodalSolutionStepVariable(Kratos.VISCOSITY)
        self.fluid_model_part.AddNodalSolutionStepVariable(
            Kratos.KINEMATIC_VISCOSITY)
        self.fluid_model_part.AddNodalSolutionStepVariable(
            Kratos.TURBULENT_VISCOSITY)
        self.fluid_model_part.AddNodalSolutionStepVariable(
            KratosRANS.RANS_Y_PLUS)

        Kratos.Logger.PrintInfo(self.__class__.__name__,
                                "Successfully added solution step variables.")

    def PrepareModelPart(self):
        self.domain_size = self.fluid_model_part.ProcessInfo[Kratos.
                                                             DOMAIN_SIZE]

        from model_part_factory import CreateDuplicateModelPart
        for element, condition in zip(self.model_elements_list,
                                      self.model_conditions_list):
            element_name = "{0}{1}D{2}N".format(element, self.domain_size,
                                                self.domain_size + 1)
            condition_name = "{0}{1}D{2}N".format(condition, self.domain_size,
                                                  self.domain_size)
            model_part = CreateDuplicateModelPart(
                self.fluid_model_part, "TurbulenceModelPart_" + element,
                element_name, condition_name)
            self.model_parts_list.append(model_part)

    def GetYPlusModel(self):
        if self.y_plus_model_process is None:
            y_plus_model_settings = self.settings["y_plus_model"]
            import rans_y_plus_model_factory
            rans_y_plus_model_factory.InitializeModelPartName(
                y_plus_model_settings, self.model, self.fluid_model_part)
            self.y_plus_model_process = rans_y_plus_model_factory.Factory(
                y_plus_model_settings, self.model)
            Kratos.Logger.PrintInfo(
                self.__class__.__name__,
                "Initialized " + self.y_plus_model_process.__str__()[:-1])

        return self.y_plus_model_process

    def GetWallVelocityModel(self):
        if self.wall_velocity_model_process is None:
            wall_velocity_model_settings = self.settings["wall_velocity_calculation"]
            import rans_wall_velocity_model_factory
            rans_wall_velocity_model_factory.InitializeModelPartName(
                wall_velocity_model_settings, self.model,
                self.fluid_model_part)
            self.wall_velocity_model_process = rans_wall_velocity_model_factory.Factory(
                wall_velocity_model_settings, self.model)
            Kratos.Logger.PrintInfo(
                self.__class__.__name__,
                "Initialized " + self.wall_velocity_model_process.__str__()[:-1])

        return self.wall_velocity_model_process

    def GetWallDistanceModel(self):
        if (self.distance_calculation_process is None):
            wall_distance_model_settings = self.settings[
                "wall_distance_calculation"]
            import rans_wall_distance_model_factory
            rans_wall_distance_model_factory.InitializeModelPartName(
                wall_distance_model_settings, self.model,
                self.fluid_model_part)
            self.distance_calculation_process = rans_wall_distance_model_factory.Factory(
                wall_distance_model_settings, self.model)

            self.distance_calculation_process.Check()
            Kratos.Logger.PrintInfo(
                self.__class__.__name__,
                "Initialized " + self.distance_calculation_process.__str__()[:-1])

        return self.distance_calculation_process

    def CreateStrategy(self, solver_settings, scheme_settings, model_part,
                       scalar_variable, scalar_variable_rate,
                       relaxed_scalar_variable_rate):
        import python_linear_solver_factory as linear_solver_factory

        default_solver_settings = Kratos.Parameters(r'''{
                "relative_tolerance"    : 1e-3,
                "absolute_tolerance"    : 1e-5,
                "max_iterations"        : 200,
                "relaxation_factor"     : 0.5,
                "echo_level"            : 0,
                "linear_solver_settings": {
                    "solver_type"  : "amgcl"
                },
                "reform_dofs_at_each_step": true,
                "move_mesh_strategy": 0,
                "move_mesh_flag": false,
                "compute_reactions": false
        }''')

        default_scheme_settings = Kratos.Parameters(r'''{
            "scheme_type": "bossak",
            "alpha_bossak": -0.3
        }''')

        solver_settings.ValidateAndAssignDefaults(default_solver_settings)
        scheme_settings.ValidateAndAssignDefaults(default_scheme_settings)

        linear_solver = linear_solver_factory.ConstructSolver(
            solver_settings["linear_solver_settings"])
        convergence_criteria = KratosRANS.GenericScalarConvergenceCriteria(
            solver_settings["relative_tolerance"].GetDouble(),
            solver_settings["absolute_tolerance"].GetDouble())
        builder_and_solver = Kratos.ResidualBasedBlockBuilderAndSolver(
            linear_solver)

        if (scheme_settings["scheme_type"].GetString() == "bossak"):
            time_scheme = KratosRANS.GenericResidualBasedBossakVelocityDynamicScalarScheme(
                scheme_settings["alpha_bossak"].GetDouble(),
                solver_settings["relaxation_factor"].GetDouble(),
                scalar_variable, scalar_variable_rate,
                relaxed_scalar_variable_rate)
        elif (scheme_settings["scheme_type"].GetString() == "steady"):
            time_scheme = KratosRANS.GenericResidualBasedSimpleSteadyScalarScheme(
                solver_settings["relaxation_factor"].GetDouble())
            self.fluid_model_part.ProcessInfo[Kratos.BOSSAK_ALPHA] = 1.0
            self.fluid_model_part.ProcessInfo[
                KratosRANS.IS_CO_SOLVING_PROCESS_ACTIVE] = True
        else:
            raise Exception("Unknown scheme_type = \"" +
                            scheme_settings["scheme_type"] + "\"")

        strategy = Kratos.ResidualBasedNewtonRaphsonStrategy(
            model_part, time_scheme, linear_solver, convergence_criteria,
            builder_and_solver, solver_settings["max_iterations"].GetInt(),
            solver_settings["compute_reactions"].GetBool(),
            solver_settings["reform_dofs_at_each_step"].GetBool(),
            solver_settings["move_mesh_flag"].GetBool())

        strategy.SetEchoLevel(solver_settings["echo_level"].GetInt() - 2)
        builder_and_solver.SetEchoLevel(
            solver_settings["echo_level"].GetInt() - 3)
        convergence_criteria.SetEchoLevel(
            solver_settings["echo_level"].GetInt() - 1)

        Kratos.Logger.PrintInfo(self.__class__.__name__,
                                "Successfully created solving strategy.")

        return strategy, linear_solver, convergence_criteria, builder_and_solver, time_scheme

    def Initialize(self):
        self.__InitializeModelPart()

        rans_variable_utils = KratosRANS.RansVariableUtils()
        rans_variable_utils.CopyScalarVar(Kratos.VISCOSITY,
                                          Kratos.KINEMATIC_VISCOSITY,
                                          self.fluid_model_part.Nodes)
        rans_variable_utils.SetScalarVar(Kratos.TURBULENT_VISCOSITY,
                                         self.nu_t_min,
                                         self.fluid_model_part.Nodes)

        self.PrepareSolvingStrategy()

        self.GetYPlusModel().ExecuteInitialize()
        self.GetWallVelocityModel().ExecuteInitialize()
        self.GetTurbulenceSolvingProcess().ExecuteInitialize()

        for strategy in self.strategies_list:
            strategy.Initialize()

        Kratos.Logger.PrintInfo(self.__class__.__name__,
                                "Initialization successfull.")

    def Check(self):
        self.GetYPlusModel().Check()
        self.GetWallVelocityModel().Check()
        self.GetTurbulenceSolvingProcess().Check()

        Kratos.Logger.PrintInfo(self.__class__.__name__,
                                "Check successfull.")

    def InitializeSolutionStep(self):
        self.GetYPlusModel().ExecuteInitializeSolutionStep()
        self.GetWallVelocityModel().ExecuteInitializeSolutionStep()
        self.GetTurbulenceSolvingProcess().ExecuteInitializeSolutionStep()

    def FinalizeSolutionStep(self):
        self.GetYPlusModel().ExecuteFinalizeSolutionStep()
        self.GetWallVelocityModel().ExecuteFinalizeSolutionStep()
        self.GetTurbulenceSolvingProcess().ExecuteFinalizeSolutionStep()

    def __InitializeModelPart(self):
        self.GetWallDistanceModel().Execute()
        Kratos.Logger.PrintInfo(self.__class__.__name__,
                                "Model part initialized.")

