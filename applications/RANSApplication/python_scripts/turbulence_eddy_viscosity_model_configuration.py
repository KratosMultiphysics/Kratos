import KratosMultiphysics as Kratos
import KratosMultiphysics.RANSApplication as KratosRANS

from KratosMultiphysics import IsDistributedRun
from KratosMultiphysics.kratos_utilities import CheckIfApplicationsAvailable
from KratosMultiphysics.process_factory import KratosProcessFactory
from KratosMultiphysics.RANSApplication.model_part_factory import CreateDuplicateModelPart

if CheckIfApplicationsAvailable("FluidDynamicsApplication"):
    import KratosMultiphysics.FluidDynamicsApplication as KratosCFD
    from KratosMultiphysics.FluidDynamicsApplication.turbulence_model_solver import TurbulenceModelSolver
else:
    msg = "RANSApplication requires FluidDynamicsApplication which is not found."
    msg += " Please install/compile it and try again."
    raise Exception(msg)

if (IsDistributedRun()
        and CheckIfApplicationsAvailable("TrilinosApplication")):
    from KratosMultiphysics.TrilinosApplication import trilinos_linear_solver_factory as linear_solver_factory
    from KratosMultiphysics.RANSApplication.TrilinosExtension import MPIGenericResidualBasedSimpleSteadyScalarScheme as steady_scheme
    from KratosMultiphysics.RANSApplication.TrilinosExtension import MPIGenericResidualBasedBossakVelocityDynamicScalarScheme as dynamic_scheme
    from KratosMultiphysics.RANSApplication.TrilinosExtension import MPIGenericScalarConvergenceCriteria as scalar_convergence_criteria
    from KratosMultiphysics.TrilinosApplication import TrilinosResidualCriteria as residual_criteria
    from KratosMultiphysics.TrilinosApplication import TrilinosNewtonRaphsonStrategy as newton_raphson_strategy
    from KratosMultiphysics.RANSApplication.block_builder_and_solvers import TrilinosPeriodicBlockBuilderAndSolver as periodic_block_builder_and_solver
    from KratosMultiphysics.RANSApplication.block_builder_and_solvers import TrilinosBlockBuilderAndSolver as block_builder_and_solver
elif (not IsDistributedRun()):
    from KratosMultiphysics import python_linear_solver_factory as linear_solver_factory
    from KratosMultiphysics.RANSApplication import GenericResidualBasedSimpleSteadyScalarScheme as steady_scheme
    from KratosMultiphysics.RANSApplication import GenericResidualBasedBossakVelocityDynamicScalarScheme as dynamic_scheme
    from KratosMultiphysics.RANSApplication import GenericScalarConvergenceCriteria as scalar_convergence_criteria
    from Kratos import ResidualCriteria as residual_criteria
    from Kratos import ResidualBasedNewtonRaphsonStrategy as newton_raphson_strategy
    from KratosMultiphysics.RANSApplication.block_builder_and_solvers import PeriodicBlockBuilderAndSolver as periodic_block_builder_and_solver
    from KratosMultiphysics.RANSApplication.block_builder_and_solvers import BlockBuilderAndSolver as block_builder_and_solver
else:
    raise Exception("Distributed run requires TrilinosApplication")


class TurbulenceEddyViscosityModelConfiguration(TurbulenceModelSolver):
    def __init__(self, model, settings):
        self._validate_settings_in_baseclass = True  # To be removed eventually

        super(TurbulenceEddyViscosityModelConfiguration, self).__init__(
            model, settings)

        # self.mesh_moving = self.settings["mesh_moving"].GetBool()
        self.strategies_list = []

        # TODO: Implement stuff for mesh_moving

        self.is_computing_solution = False
        self.EpetraCommunicator = None

        self.model_elements_list = []
        self.model_conditions_list = []
        self.model_parts_list = []

    def SetCommunicator(self, communicator):
        self.EpetraCommunicator = communicator

    def GetDefaultSettings(self):
        return Kratos.Parameters(r'''{
            "model_type"            : "",
            "velocity_pressure_relaxation_factor": 0.2,
            "model_settings"        : {},
            "auxiliar_process_list"   : [],
            "mesh_moving"             : false,
            "echo_level"              : 0
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

    def SetParentSolvingStrategy(self, parent_solving_strategy):
        self.GetTurbulenceSolvingProcess().SetParentSolvingStrategy(parent_solving_strategy)

    def __CreateBuilderAndSolver(self, linear_solver, is_periodic):
        if (is_periodic):
            return periodic_block_builder_and_solver(linear_solver, self.EpetraCommunicator)
        else:
            return block_builder_and_solver(linear_solver,  self.EpetraCommunicator)

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

        linear_solver = linear_solver_factory.ConstructSolver(
            solver_settings["linear_solver_settings"])

        is_periodic = solver_settings["is_periodic"].GetBool()

        if is_periodic:
            self.__InitializePeriodicConditions(model_part, scalar_variable)

        # TODO:
        if is_periodic and IsDistributedRun():
            msg = "\nCurrently periodic conditions in mpi is not supported due to following reasons:\n\n"
            msg += "    1. TrilinosResidualCriteria [ConvergenceCriterian]\n"
            msg += "PeriodicConditions duplicates one patch's equation ids to the counter patch. "
            msg += "The node and its corresponding dof might not fall in to the same partition raising an error in convergence calculation.\n\n"
            msg += "    2. ConnectivityPreserveModeller\n"
            msg += "Currently connectivity preserve modeller replaces all the conditions in an mdpa with given new condition. "
            msg += "This modeller is used to create modelparts having k-epsilon elements and conditions while sharing the same nodes as in VMS solution. "
            msg += "In the case of MPI, it is essential to have the PeriodicConditions in the mdpa file in order to properly distribute nodes to partitions using MetisApplication. "
            msg += "But if this is the case, PeriodicConditions also will be replaced by k-epsilon specific conditions casuing a segmentation fault.\n"
            msg += "    3. TrilinosBlockBuilderAndSolverPeriodic\n"
            msg += "In the case of MPI periodic in 2D, problem uses TrilinosBlockBuilderAndSolverPeriodic block builder and solver, which identifies "
            msg += "periodic conditions by number of nodes in the condition. So, In 2D all wall conditions and PeriodicConditions have only 2 nodes, all will be "
            msg += "considered as PeriodicConditions and will make the global assembly accordingly which is wrong."
            msg += "Therefore this error msg is printed in order to avoid confusion."
            raise Exception(msg)

        builder_and_solver = self.__CreateBuilderAndSolver(
            linear_solver, is_periodic)

        if (scheme_settings["scheme_type"].GetString() == "bossak"):
            convergence_criteria_type = scalar_convergence_criteria
        elif (scheme_settings["scheme_type"].GetString() == "steady"):
            convergence_criteria_type = residual_criteria

        convergence_criteria = convergence_criteria_type(
            solver_settings["relative_tolerance"].GetDouble(),
            solver_settings["absolute_tolerance"].GetDouble())

        if (scheme_settings["scheme_type"].GetString() == "bossak"):
            time_scheme = dynamic_scheme(
                scheme_settings["alpha_bossak"].GetDouble(),
                solver_settings["relaxation_factor"].GetDouble(),
                scalar_variable, scalar_variable_rate,
                relaxed_scalar_variable_rate)
        elif (scheme_settings["scheme_type"].GetString() == "steady"):
            time_scheme = steady_scheme(
                solver_settings["relaxation_factor"].GetDouble())
            self.fluid_model_part.ProcessInfo[
                Kratos.
                BOSSAK_ALPHA] = scheme_settings["alpha_bossak"].GetDouble()
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

        strategy = newton_raphson_strategy(
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

        return strategy

    def Initialize(self):
        self.PrepareSolvingStrategy()

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

    def Finalize(self):
        for strategy in self.strategies_list:
            strategy.Clear()

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
