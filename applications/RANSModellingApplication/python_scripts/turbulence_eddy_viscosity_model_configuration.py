import KratosMultiphysics as Kratos
import KratosMultiphysics.RANSModellingApplication as KratosRANS
import math

from KratosMultiphysics.kratos_utilities import CheckIfApplicationsAvailable
from KratosMultiphysics.RANSModellingApplication.model_part_factory import CreateDuplicateModelPart

if CheckIfApplicationsAvailable("FluidDynamicsApplication"):
    import KratosMultiphysics.FluidDynamicsApplication as KratosCFD
    from KratosMultiphysics.FluidDynamicsApplication.turbulence_model_configuration import TurbulenceModelConfiguration
else:
    msg = "RANSModellingApplication requires FluidDynamicsApplication which is not found."
    msg += " Please install/compile it and try again."
    raise Exception(msg)

if CheckIfApplicationsAvailable("TrilinosApplication"):
    import KratosMultiphysics.mpi as KratosMPI  # MPI-python interface
    # Import applications
    import KratosMultiphysics.TrilinosApplication as KratosTrilinos  # MPI solvers

class TurbulenceEddyViscosityModelConfiguration(TurbulenceModelConfiguration):
    def __init__(self, model, settings):
        self._validate_settings_in_baseclass = True  # To be removed eventually
        self.parallel_type = ""

        super(TurbulenceEddyViscosityModelConfiguration, self).__init__(
            model, settings)

        # self.mesh_moving = self.settings["mesh_moving"].GetBool()
        self.turbulence_model_process = None
        self.strategies_list = []
        self.nu_t_min = self.settings["turbulent_viscosity_min"].GetDouble()
        self.nu_t_max = self.settings["turbulent_viscosity_max"].GetDouble()

        # TODO: Implement stuff for mesh_moving

        self.is_computing_solution = False

        self.model_elements_list = []
        self.model_conditions_list = []
        self.model_parts_list = []

    def SetParallelType(self, parallel_type):
        self.parallel_type = parallel_type

        if (self.parallel_type == "MPI"):
            if CheckIfApplicationsAvailable("TrilinosApplication"):
                import KratosMultiphysics.mpi as KratosMPI  # MPI-python interface
                # Import applications
                import KratosMultiphysics.TrilinosApplication as KratosTrilinos  # MPI solvers
            else:
                raise Exception(
                    "MPI parallel type enforced without TrilinosApplication. Please install/compile it and try again."
                )

    def SetCommunicator(self, communicator):
            self.EpetraCommunicator = communicator

    def GetDefaultSettings(self):
        return Kratos.Parameters(r'''{
            "model_type"            : "",
            "velocity_pressure_relaxation_factor": 0.2,
            "model_settings"        : {},
            "auxiliar_process_list"   : [],
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
        self.fluid_model_part.AddNodalSolutionStepVariable(
            Kratos.RELAXED_ACCELERATION)

        Kratos.Logger.PrintInfo(self.__class__.__name__,
                                "Successfully added solution step variables.")

    def PrepareModelPart(self):
        self.domain_size = self.fluid_model_part.ProcessInfo[Kratos.
                                                             DOMAIN_SIZE]

        original_condition_name = self.GetFluidVelocityPressureConditionName()

        for element, condition in zip(self.model_elements_list,
                                      self.model_conditions_list):
            element_name = "{0}{1}D{2}N".format(element, self.domain_size,
                                                self.domain_size + 1)
            condition_name = "{0}{1}D{2}N".format(condition, self.domain_size,
                                                  self.domain_size)
            model_part = CreateDuplicateModelPart(
                self.fluid_model_part, "TurbulenceModelPart_" + element,
                element_name, condition_name, original_condition_name)
            self.model_parts_list.append(model_part)

    def __CreateLinearSolver(self, linear_solver_settings):
        if (self.parallel_type == "MPI"):
            from KratosMultiphysics.TrilinosApplication import trilinos_linear_solver_factory as linear_solver_factory
        elif (self.parallel_type == "OpenMP"):
            from KratosMultiphysics import python_linear_solver_factory as linear_solver_factory

        return linear_solver_factory.ConstructSolver(linear_solver_settings)

    def __CreateBuilderAndSolver(self, linear_solver, is_periodic):
        if (self.parallel_type == "MPI"):
            if (is_periodic):
                return KratosTrilinos.TrilinosBlockBuilderAndSolverPeriodic(
                    self.EpetraCommunicator, 30,
                    linear_solver, KratosCFD.PATCH_INDEX)
            else:
                return KratosTrilinos.TrilinosBlockBuilderAndSolver(
                    self.EpetraCommunicator, 30,
                    linear_solver)
        elif (self.parallel_type == "OpenMP"):
            if (is_periodic):
                return KratosCFD.ResidualBasedBlockBuilderAndSolverPeriodic(
                    linear_solver, KratosCFD.PATCH_INDEX)
            else:
                return Kratos.ResidualBasedBlockBuilderAndSolver(linear_solver)

    def __CreateConvergenceCriteria(self, scheme_type, relative_tolerance,
                                    absolute_tolerance, is_periodic):
        if (self.parallel_type == "MPI"):
            if (scheme_type == "bossak"):
                return KratosRANS.MPIGenericScalarConvergenceCriteria(
                    relative_tolerance, absolute_tolerance)
            elif (scheme_type == "steady"):
                return KratosTrilinos.TrilinosResidualCriteria(relative_tolerance,
                                                       absolute_tolerance)
        elif (self.parallel_type == "OpenMP"):
            if (scheme_type == "bossak"):
                return KratosRANS.GenericScalarConvergenceCriteria(
                    relative_tolerance, absolute_tolerance)
            elif (scheme_type == "steady"):
                return Kratos.ResidualCriteria(relative_tolerance,
                                               absolute_tolerance)

    def __CreateDynamicTimeScheme(self, alpha_bossak, relaxation_factor,
                                  scalar_variable, scalar_variable_rate,
                                  relaxed_scalar_variable_rate):
        if (self.parallel_type == "MPI"):
            return KratosRANS.MPIGenericResidualBasedBossakVelocityDynamicScalarScheme(
                alpha_bossak, relaxation_factor, scalar_variable,
                scalar_variable_rate, relaxed_scalar_variable_rate)
        elif (self.parallel_type == "OpenMP"):
            return KratosRANS.GenericResidualBasedBossakVelocityDynamicScalarScheme(
                alpha_bossak, relaxation_factor, scalar_variable,
                scalar_variable_rate, relaxed_scalar_variable_rate)

    def __CreateSteadyScheme(self, relaxation_factor):
        if (self.parallel_type == "MPI"):
            return KratosRANS.MPIGenericResidualBasedSimpleSteadyScalarScheme(
                relaxation_factor)
        elif (self.parallel_type == "OpenMP"):
            return KratosRANS.GenericResidualBasedSimpleSteadyScalarScheme(
                relaxation_factor)

    def CreateStrategy(self, solver_settings, scheme_settings, model_part,
                       scalar_variable, scalar_variable_rate,
                       relaxed_scalar_variable_rate):
        default_solver_settings = Kratos.Parameters(r'''{
                "is_periodic"           : false,
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

        linear_solver = self.__CreateLinearSolver(
            solver_settings["linear_solver_settings"])

        is_periodic = solver_settings["is_periodic"].GetBool()

        if is_periodic:
            self.__InitializePeriodicConditions(model_part, scalar_variable)

        builder_and_solver = self.__CreateBuilderAndSolver(
            linear_solver, is_periodic)

        convergence_criteria = self.__CreateConvergenceCriteria(
            scheme_settings["scheme_type"].GetString(),
            solver_settings["relative_tolerance"].GetDouble(),
            solver_settings["absolute_tolerance"].GetDouble(), is_periodic)

        if (scheme_settings["scheme_type"].GetString() == "bossak"):
            time_scheme = self.__CreateDynamicTimeScheme(
                scheme_settings["alpha_bossak"].GetDouble(),
                solver_settings["relaxation_factor"].GetDouble(),
                scalar_variable, scalar_variable_rate,
                relaxed_scalar_variable_rate)
        elif (scheme_settings["scheme_type"].GetString() == "steady"):
            time_scheme = self.__CreateSteadyScheme(
                solver_settings["relaxation_factor"].GetDouble())
            self.fluid_model_part.ProcessInfo[
                Kratos.
                BOSSAK_ALPHA] = scheme_settings["alpha_bossak"].GetDouble()
            self.fluid_model_part.ProcessInfo[Kratos.BOSSAK_ALPHA] = 1.0
            self.fluid_model_part.ProcessInfo[
                KratosRANS.IS_CO_SOLVING_PROCESS_ACTIVE] = True
            if (self.fluid_model_part.ProcessInfo[Kratos.DYNAMIC_TAU] != 0.0):
                Kratos.Logger.PrintWarning(
                    self.__class__.__name__,
                    "Steady solution doesn't have zero DYNAMIC_TAU [ DYNAMIC_TAU = "
                    + str(
                        self.fluid_model_part.ProcessInfo[Kratos.DYNAMIC_TAU])
                    + " ].")
        else:
            raise Exception("Unknown scheme_type = \"" +
                            scheme_settings["scheme_type"] + "\"")

        if (self.parallel_type == "MPI"):
            strategy_type = KratosTrilinos.TrilinosNewtonRaphsonStrategy
        else:
            strategy_type = Kratos.ResidualBasedNewtonRaphsonStrategy

        strategy = strategy_type(
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

        if (is_periodic):
            Kratos.Logger.PrintInfo(
                self.__class__.__name__,
                "Successfully created periodic solving strategy for " +
                scalar_variable.Name() + ".")
        else:
            Kratos.Logger.PrintInfo(
                self.__class__.__name__,
                "Successfully created solving strategy for " +
                scalar_variable.Name() + ".")

        return strategy, linear_solver, convergence_criteria, builder_and_solver, time_scheme

    def Initialize(self):
        rans_variable_utils = KratosRANS.RansVariableUtils()
        # set this from the json file.
        # rans_variable_utils.CopyScalarVar(Kratos.VISCOSITY,
        #                                   Kratos.KINEMATIC_VISCOSITY,
        #                                   self.fluid_model_part.Nodes)
        # rans_variable_utils.SetScalarVar(Kratos.TURBULENT_VISCOSITY,
        #                                  self.nu_t_min,
        #                                  self.fluid_model_part.Nodes)

        self.PrepareSolvingStrategy()

        from process_factory import KratosProcessFactory
        factory = KratosProcessFactory(self.model)
        self.auxiliar_process_list = factory.ConstructListOfProcesses(
            self.settings["auxiliar_process_list"])
        for process in self.auxiliar_process_list:
            self.GetTurbulenceSolvingProcess().AddAuxiliaryProcess(process)
            process.ExecuteInitialize()

        self.GetTurbulenceSolvingProcess().ExecuteInitialize()

        for strategy in self.strategies_list:
            strategy.Initialize()

        Kratos.Logger.PrintInfo(self.__class__.__name__,
                                "Initialization successfull.")

    def Check(self):
        self.GetTurbulenceSolvingProcess().Check()

        for process in self.auxiliar_process_list:
            process.Check()

        Kratos.Logger.PrintInfo(self.__class__.__name__, "Check successfull.")

    def InitializeSolutionStep(self):
        for process in self.auxiliar_process_list:
            process.ExecuteInitializeSolutionStep()

        self.GetTurbulenceSolvingProcess().ExecuteInitializeSolutionStep()

    def FinalizeSolutionStep(self):
        for process in self.auxiliar_process_list:
            process.ExecuteFinalizeSolutionStep()

        self.GetTurbulenceSolvingProcess().ExecuteFinalizeSolutionStep()

    def __InitializePeriodicConditions(self, model_part, scalar_variable):
        properties = model_part.CreateNewProperties(
            model_part.NumberOfProperties() + 1)
        pcu = KratosCFD.PeriodicConditionUtilities(
            model_part, model_part.ProcessInfo[Kratos.DOMAIN_SIZE])
        pcu.AddPeriodicVariable(properties, scalar_variable)

        index = model_part.NumberOfConditions()
        for condition in self.fluid_model_part.Conditions:
            if condition.Is(Kratos.PERIODIC):
                index += 1
                node_id_list = [node.Id for node in condition.GetNodes()]
                periodic_condition = model_part.CreateNewCondition(
                    "PeriodicCondition", index, node_id_list, properties)
                periodic_condition.Set(Kratos.PERIODIC)
