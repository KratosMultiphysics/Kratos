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
        default_settings = Kratos.Parameters(r'''{
            "model_type"            : "",
            "model_settings"        : {},
            "distance_calculation"  : {
                "max_iterations"         : 5,
                "linear_solver_settings" : {}
            },
            "y_plus_model"      : {
                "model_type"     : "",
                "model_settings" : {}
            },
            "mesh_moving"             : false,
            "echo_level"              : 0,
            "turbulent_viscosity_min" : 1e-12,
            "turbulent_viscosity_max" : 1e+2
        }''')
        self.settings = settings.ValidateAndAssignDefaults(default_settings)

        super(TurbulenceEddyViscosityModelConfiguration, self).__init__(model, settings)

        # self.mesh_moving = self.settings["mesh_moving"].GetBool()
        self.distance_calculation_process = None
        self.y_plus_model_process = None
        self.turbulence_model_process = None
        self.strategies_list = []
        self.nu_t_min = self.settings["turbulent_viscosity_min"].GetDouble()
        self.nu_t_max = self.settings["turbulent_viscosity_max"].GetDouble()

        # TODO: Implement stuff for mesh_moving

        self.is_computing_solution = False

        self.model_elements_list = []
        self.model_conditions_list = []
        self.model_parts_list = []

    def Initialize(self):
        self.__InitializeModelPart()

        super(TurbulenceEddyViscosityModelConfiguration, self).Initialize()

        rans_variable_utils = KratosRANS.RansVariableUtils()
        rans_variable_utils.CopyScalarVar(Kratos.VISCOSITY, Kratos.KINEMATIC_VISCOSITY,  self.fluid_model_part.Nodes)
        rans_variable_utils.SetScalarVar(Kratos.TURBULENT_VISCOSITY, self.nu_t_min, self.fluid_model_part.Nodes)

        self.PrepareSolvingStrategy()

        for strategy in self.strategies_list:
            strategy.Initialize()

        Kratos.Logger.PrintInfo(self.__class__.__name__, "Initialization successfull.")

    def AddVariables(self):
        # adding variables required by rans eddy viscosity models
        self.fluid_model_part.AddNodalSolutionStepVariable(Kratos.DISTANCE)
        self.fluid_model_part.AddNodalSolutionStepVariable(Kratos.FLAG_VARIABLE)
        self.fluid_model_part.AddNodalSolutionStepVariable(Kratos.VISCOSITY)
        self.fluid_model_part.AddNodalSolutionStepVariable(Kratos.KINEMATIC_VISCOSITY)
        self.fluid_model_part.AddNodalSolutionStepVariable(Kratos.TURBULENT_VISCOSITY)
        self.fluid_model_part.AddNodalSolutionStepVariable(KratosRANS.RANS_Y_PLUS)
        self.fluid_model_part.AddNodalSolutionStepVariable(KratosRANS.OLD_CONVERGENCE_VARIABLE)

        Kratos.Logger.PrintInfo(self.__class__.__name__, "Successfully added solution step variables.")

    def Check(self):
        self.GetYPlusModel().Check()
        self.GetTurbulenceSolvingProcess().Check()

    def PrepareModelPart(self):
        self.domain_size = self.fluid_model_part.ProcessInfo[Kratos.DOMAIN_SIZE]

        from model_part_factory import CreateDuplicateModelPart
        for element, condition in zip(self.model_elements_list, self.model_conditions_list):
            element_name = "{0}{1}D{2}N".format(element, self.domain_size,
                                                self.domain_size + 1)
            condition_name = "{0}{1}D{2}N".format(
                condition, self.domain_size, self.domain_size)
            model_part = CreateDuplicateModelPart(self.fluid_model_part, "TurbulenceModelPart_" + element, element_name, condition_name)
            self.model_parts_list.append(model_part)


    def GetYPlusModel(self):
        if self.y_plus_model_process is None:
            import rans_y_plus_model_factory
            self.y_plus_model_process = rans_y_plus_model_factory.Factory(
                self.fluid_model_part, self.settings["y_plus_model"])
            Kratos.Logger.PrintInfo(self.__class__.__name__,
                      "Initialized " + self.y_plus_model_process.__str__())

        return self.y_plus_model_process

    def CreateStrategy(self, solver_settings, scheme_settings, model_part,
                       scalar_variable, scalar_variable_rate,
                       relaxed_scalar_variable_rate):
        import python_linear_solver_factory as linear_solver_factory
        linear_solver = linear_solver_factory.ConstructSolver(
            solver_settings["linear_solver_settings"])
        convergence_criteria = KratosRANS.GenericScalarConvergenceCriteria(
            solver_settings["relative_tolerance"].GetDouble(),
            solver_settings["absolute_tolerance"].GetDouble())
        builder_and_solver = Kratos.ResidualBasedBlockBuilderAndSolver(
            linear_solver)
        time_scheme = KratosRANS.GenericResidualBasedBossakVelocityDynamicScalarScheme(
            scheme_settings["alpha_bossak"].GetDouble(), scalar_variable,
            scalar_variable_rate, relaxed_scalar_variable_rate)

        strategy = Kratos.ResidualBasedNewtonRaphsonStrategy(
            model_part, time_scheme, linear_solver, convergence_criteria,
            builder_and_solver, solver_settings["max_iterations"].GetInt(),
            False, False, False)

        strategy.SetEchoLevel(solver_settings["echo_level"].GetInt())
        builder_and_solver.SetEchoLevel(solver_settings["echo_level"].GetInt())
        convergence_criteria.SetEchoLevel(solver_settings["echo_level"].GetInt())

        Kratos.Logger.PrintInfo(self.__class__.__name__, "Successfully created solving strategy.")

        return strategy, linear_solver, convergence_criteria, builder_and_solver, time_scheme

    def __InitializeModelPart(self):
        variable_utils = Kratos.VariableUtils()
        variable_utils.SetScalarVar(Kratos.DISTANCE, 1.0, self.fluid_model_part.Nodes)
        variable_utils.SetScalarVar(Kratos.DISTANCE, 0.0, self.fluid_model_part.Nodes, Kratos.STRUCTURE, True)

        self.__CalculateWallDistances()

        Kratos.Logger.PrintInfo(self.__class__.__name__, "Model part initialized.")

    def __CalculateWallDistances(self):
        if (self.distance_calculation_process is None):
            import python_linear_solver_factory as linear_solver_factory
            self.distance_calculation_linear_solver = linear_solver_factory.ConstructSolver(
                self.settings["distance_calculation"]
                ["linear_solver_settings"])
            max_iterations = self.settings["distance_calculation"][
                "max_iterations"].GetInt()
            if (self.domain_size == 2):
                self.distance_calculation_process = Kratos.VariationalDistanceCalculationProcess2D(
                    self.fluid_model_part,
                    self.distance_calculation_linear_solver, max_iterations)
            elif (self.domain_size == 3):
                self.distance_calculation_process = Kratos.VariationalDistanceCalculationProcess3D(
                    self.fluid_model_part,
                    self.distance_calculation_linear_solver, max_iterations)
            else:
                raise Exception("Unsupported domain size")

            Kratos.Logger.PrintInfo(
                self.__class__.__name__,
                "Variational distance calculation process initialized.")

        self.distance_calculation_process.Execute()
        Kratos.Logger.PrintInfo(self.__class__.__name__, "Wall distances calculated.")
