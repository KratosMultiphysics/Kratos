from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# importing the Kratos Library
from KratosMultiphysics import *
from KratosMultiphysics.FluidDynamicsApplication import *
# Check that KratosMultiphysics was imported in the main script
CheckForPreviousImport()


def AddVariables(model_part, config=None):
    model_part.AddNodalSolutionStepVariable(VELOCITY)
    model_part.AddNodalSolutionStepVariable(ACCELERATION)
    model_part.AddNodalSolutionStepVariable(MESH_VELOCITY)
    model_part.AddNodalSolutionStepVariable(PRESSURE)
    model_part.AddNodalSolutionStepVariable(IS_STRUCTURE)
    model_part.AddNodalSolutionStepVariable(DISPLACEMENT)
    model_part.AddNodalSolutionStepVariable(VISCOSITY)
    model_part.AddNodalSolutionStepVariable(DENSITY)
    model_part.AddNodalSolutionStepVariable(BODY_FORCE)
    model_part.AddNodalSolutionStepVariable(NODAL_AREA)
    model_part.AddNodalSolutionStepVariable(NODAL_H)
    model_part.AddNodalSolutionStepVariable(ADVPROJ)
    model_part.AddNodalSolutionStepVariable(DIVPROJ)
    model_part.AddNodalSolutionStepVariable(REACTION)
    model_part.AddNodalSolutionStepVariable(REACTION_WATER_PRESSURE)
    model_part.AddNodalSolutionStepVariable(EXTERNAL_PRESSURE)
    model_part.AddNodalSolutionStepVariable(FLAG_VARIABLE)
    model_part.AddNodalSolutionStepVariable(NORMAL)
    model_part.AddNodalSolutionStepVariable(Y_WALL)
    model_part.AddNodalSolutionStepVariable(PATCH_INDEX)

    if config is not None:
        if hasattr(config, "TurbulenceModel"):
            if config.TurbulenceModel == "Spalart-Allmaras":
                model_part.AddNodalSolutionStepVariable(TURBULENT_VISCOSITY)
                model_part.AddNodalSolutionStepVariable(MOLECULAR_VISCOSITY)
                model_part.AddNodalSolutionStepVariable(TEMP_CONV_PROJ)
                model_part.AddNodalSolutionStepVariable(DISTANCE)

    print("variables for the vms monolithic solver added correctly")


def AddDofs(model_part, config=None):
    for node in model_part.Nodes:
        # adding dofs
        node.AddDof(VELOCITY_X, REACTION_X)
        node.AddDof(VELOCITY_Y, REACTION_Y)
        node.AddDof(VELOCITY_Z, REACTION_Z)
        node.AddDof(PRESSURE, REACTION_WATER_PRESSURE)

    if config is not None:
        if hasattr(config, "TurbulenceModel"):
            if config.TurbulenceModel == "Spalart-Allmaras":
                for node in model_part.Nodes:
                    node.AddDof(TURBULENT_VISCOSITY)

    print("dofs for the vms monolithic solver added correctly")


class MonolithicSolver:

    def __setattr__(self, name, value):
        if name == "periodic":
            raise TypeError("'periodic' attribute no longer supported.")
        else:
            self.__dict__[name] = value

    def __init__(self, model_part, domain_size):

        self.model_part = model_part
        self.domain_size = domain_size

        self.alpha = -0.3
        self.move_mesh_strategy = 0

        # definition of the solvers
        try:
            from KratosMultiphysics.ExternalSolversApplication import SuperLUIterativeSolver
            self.linear_solver = SuperLUIterativeSolver()
        except:
            self.linear_solver = SkylineLUFactorizationSolver()

        # definition of the convergence criteria
        self.rel_vel_tol = 1e-5
        self.abs_vel_tol = 1e-7
        self.rel_pres_tol = 1e-5
        self.abs_pres_tol = 1e-7

        self.dynamic_tau = 0.0
        self.oss_switch = 0

        # non newtonian setting
        self.regularization_coef = 1000

        self.max_iter = 30

        # default settings
        self.echo_level = 0
        self.compute_reactions = True
        self.ReformDofSetAtEachStep = True
        self.CalculateNormDxFlag = True
        self.MoveMeshFlag = False
        self.use_slip_conditions = False
        self.time_scheme = None
        self.builder_and_solver = None
        self.turbulence_model = None
        self.use_spalart_allmaras = False
        self.use_des = False
        self.Cdes = 1.0
        self.wall_nodes = list()
        self.spalart_allmaras_linear_solver = None

        self.divergence_clearance_steps = 0

        print("Construction monolithic solver finished")

    #
    def Initialize(self):
        # check if slip conditions are defined
        if self.use_slip_conditions == False:
            for cond in self.model_part.Conditions:
                if cond.GetValue(IS_STRUCTURE) != 0.0:
                    self.use_slip_conditions = True
                    break

        # if we use slip conditions, calculate normals on the boundary
        if self.use_slip_conditions:
            self.normal_util = NormalCalculationUtils()
            self.normal_util.CalculateOnSimplex(
                self.model_part, self.domain_size, IS_STRUCTURE)

            for cond in self.model_part.Conditions:
                if cond.GetValue(IS_STRUCTURE) != 0.0:
                    for node in cond.GetNodes():
                        node.SetValue(IS_STRUCTURE, 1.0)

        # Turbulence model
        if self.use_spalart_allmaras:
            self.activate_spalart_allmaras()

        # creating the solution strategy
        self.conv_criteria = VelPrCriteria(self.rel_vel_tol, self.abs_vel_tol,
                                           self.rel_pres_tol, self.abs_pres_tol)

        (self.conv_criteria).SetEchoLevel(self.echo_level)

        # custom schemes should be defined in the calling python script
        # before calling Initialize().
        #
        # Examples:
        # fluid_solver.time_scheme = ResidualBasedPredictorCorrectorVelocityBossakSchemeTurbulent(
        #                             fluid_solver.alpha,
        #                             fluid_solver.move_mesh_strategy,
        #                             domain_size,
        #                             PATCH_INDEX)
        #
        # fluid_solver.activate_spalart_allmaras()
        # fluid_solver.time_scheme = ResidualBasedSimpleSteadyScheme(
        #                             0.1, 0.1,
        #                             domain_size,
        #                             fluid_solver.turbulence_model)
        #
        # a scheme_factory may be added in the future to parse configuration
        # parameters set with a GiD problem type. below is the default scheme.
        if self.time_scheme is None:
            if self.turbulence_model is None:
                self.time_scheme = ResidualBasedPredictorCorrectorVelocityBossakSchemeTurbulent(self.alpha,
                                                                                                self.move_mesh_strategy,
                                                                                                self.domain_size)
            else:
                self.time_scheme = ResidualBasedPredictorCorrectorVelocityBossakSchemeTurbulent(self.alpha,
                                                                                                self.move_mesh_strategy,
                                                                                                self.domain_size,
                                                                                                self.turbulence_model)

        # custom builder and solvers should be defined before calling
        # Initialize().
        #
        # Example:
        # fluid_solver.builder_and_solver = ResidualBasedBlockBuilderAndSolverPeriodic(
        #                                    fluid_solver.linear_solver,
        #                                    PATCH_INDEX)
        if self.builder_and_solver is None:
            self.builder_and_solver = ResidualBasedBlockBuilderAndSolver(self.linear_solver)

        self.solver = ResidualBasedNewtonRaphsonStrategy(
            self.model_part, self.time_scheme, self.linear_solver, self.conv_criteria,
            self.builder_and_solver, self.max_iter, self.compute_reactions, self.ReformDofSetAtEachStep, self.MoveMeshFlag)
        (self.solver).SetEchoLevel(self.echo_level)
        self.solver.Check()

        self.model_part.ProcessInfo.SetValue(DYNAMIC_TAU, self.dynamic_tau)
        self.model_part.ProcessInfo.SetValue(OSS_SWITCH, self.oss_switch)
        self.model_part.ProcessInfo.SetValue(M, self.regularization_coef)

        print ("Initialization monolithic solver finished")

    def Solve(self):

        if self.divergence_clearance_steps > 0:
            print("Calculating divergence-free initial condition")
            # initialize with a Stokes solution step
            try:
                import KratosMultiphysics.ExternalSolversApplication as kes
                smoother_type = kes.AMGCLSmoother.DAMPED_JACOBI
                solver_type = kes.AMGCLIterativeSolverType.CG
                gmres_size = 50
                max_iter = 200
                tol = 1e-7
                verbosity = 0
                stokes_linear_solver = kes.AMGCLSolver(
                    smoother_type,
                    solver_type,
                    tol,
                    max_iter,
                    verbosity,
                    gmres_size)
            except:
                pPrecond = DiagonalPreconditioner()
                stokes_linear_solver = BICGSTABSolver(1e-9, 5000, pPrecond)
            stokes_process = StokesInitializationProcess(
                self.model_part,
                stokes_linear_solver,
                self.domain_size,
                PATCH_INDEX)
            # copy periodic conditions to Stokes problem
            stokes_process.SetConditions(self.model_part.Conditions)
            # execute Stokes process
            stokes_process.Execute()
            stokes_process = None

            for node in self.model_part.Nodes:
                node.SetSolutionStepValue(PRESSURE, 0, 0.0)
                node.SetSolutionStepValue(ACCELERATION_X, 0, 0.0)
                node.SetSolutionStepValue(ACCELERATION_Y, 0, 0.0)
                node.SetSolutionStepValue(ACCELERATION_Z, 0, 0.0)
#                vel = node.GetSolutionStepValue(VELOCITY)
#                for i in range(0,2):
#                    node.SetSolutionStepValue(VELOCITY,i,vel)

            self.divergence_clearance_steps = 0
            print("Finished divergence clearance")

        if self.ReformDofSetAtEachStep:
            if self.use_slip_conditions:
                self.normal_util.CalculateOnSimplex(
                    self.model_part, self.domain_size, IS_STRUCTURE)
            if self.use_spalart_allmaras:
                self.neighbour_search.Execute()

        (self.solver).Solve()

    #
    def SetEchoLevel(self, level):
        (self.solver).SetEchoLevel(level)

    #
    def Clear(self):
        (self.solver).Clear()

    #
    def activate_smagorinsky(self, C):
        for elem in self.model_part.Elements:
            elem.SetValue(C_SMAGORINSKY, C)

    #
    def activate_spalart_allmaras(self):

        # Spalart-Allmaras initialization
        for node in self.wall_nodes:
            node.SetValue(IS_VISITED, 1.0)
            node.SetSolutionStepValue(DISTANCE, 0, 0.0)

        if(self.domain_size == 2):
            self.redistance_utils = ParallelDistanceCalculator2D()
        else:
            self.redistance_utils = ParallelDistanceCalculator3D()

        max_levels = 100
        max_distance = 1000
        self.redistance_utils.CalculateDistancesLagrangianSurface(
            self.model_part,
            DISTANCE,
            NODAL_AREA,
            max_levels,
            max_distance)

        # Spalart-Allmaras uses the componentwise builder and solver, which
        # requires a neighbour search
        number_of_avg_elems = 10
        number_of_avg_nodes = 10
        self.neighbour_search = FindNodalNeighboursProcess(
            self.model_part, number_of_avg_elems, number_of_avg_nodes)
        self.neighbour_search.Execute()

        non_linear_tol = 0.001
        max_it = 50 #10
        reform_dofset = self.ReformDofSetAtEachStep
        time_order = 2

        if self.spalart_allmaras_linear_solver is None:
            #pPrecond = DiagonalPreconditioner()
            #self.spalart_allmaras_linear_solver = BICGSTABSolver(
                #1e-9, 300, pPrecond)

            import KratosMultiphysics.ExternalSolversApplication as kes
            smoother_type = kes.AMGCLSmoother.ILU0
            solver_type = kes.AMGCLIterativeSolverType.GMRES
            gmres_size = 200
            max_iter = 200
            tol = 1e-9
            verbosity = 1
            linear_solver = kes.AMGCLSolver(
                    smoother_type,
                    solver_type,
                    tol,
                    max_iter,
                    verbosity,
                    gmres_size)
            #self.spalart_allmaras_linear_solver =  ScalingSolver(linear_solver, True)

            #self.spalart_allmaras_linear_solver = kes.SuperLUSolver()
            self.spalart_allmaras_linear_solver = kes.PastixSolver(0,False)

        self.turbulence_model = SpalartAllmarasTurbulenceModel(
            self.model_part,
            self.spalart_allmaras_linear_solver,
            self.domain_size,
            non_linear_tol,
            max_it,
            reform_dofset,
            time_order)
        if self.use_des:
            self.turbulence_model.ActivateDES(self.Cdes)

#
#
def CreateSolver(model_part, config):
    fluid_solver = MonolithicSolver(model_part, config.domain_size)

    if(hasattr(config, "alpha")):
        fluid_solver.alpha = config.alpha

    # definition of the convergence criteria
    if(hasattr(config, "velocity_relative_tolerance")):
        fluid_solver.rel_vel_tol = config.velocity_relative_tolerance
    if(hasattr(config, "velocity_absolute_tolerance")):
        fluid_solver.abs_vel_tol = config.velocity_absolute_tolerance
    if(hasattr(config, "pressure_relative_tolerance")):
        fluid_solver.rel_pres_tol = config.pressure_relative_tolerance
    if(hasattr(config, "pressure_absolute_tolerance")):
        fluid_solver.abs_pres_tol = config.pressure_absolute_tolerance
    if(hasattr(config, "dynamic_tau")):
        fluid_solver.dynamic_tau = config.dynamic_tau
    if(hasattr(config, "oss_switch")):
        fluid_solver.oss_switch = config.oss_switch
    if(hasattr(config, "max_iteration")):
        fluid_solver.max_iter = config.max_iteration
    if(hasattr(config, "echo_level")):
        fluid_solver.echo_level = config.echo_level
    if(hasattr(config, "compute_reactions")):
        fluid_solver.compute_reactions = config.compute_reactions
    if(hasattr(config, "ReformDofSetAtEachStep")):
        fluid_solver.ReformDofSetAtEachStep = config.ReformDofSetAtEachStep
    if(hasattr(config, "divergence_cleareance_step")):
        fluid_solver.divergence_clearance_steps = config.divergence_cleareance_step
    if hasattr(config, "TurbulenceModel"):
        if config.TurbulenceModel == "Spalart-Allmaras":
            fluid_solver.use_spalart_allmaras = True
        elif config.TurbulenceModel == "Smagorinsky-Lilly":
            if hasattr(config, "SmagorinskyConstant"):
                fluid_solver.activate_smagorinsky(config.SmagorinskyConstant)
            else:
                msg = """Fluid solver error: Smagorinsky model requested, but
                         the value for the Smagorinsky constant is
                         undefined."""
                raise Exception(msg)

    import linear_solver_factory
    if(hasattr(config, "linear_solver_config")):
        fluid_solver.linear_solver = linear_solver_factory.ConstructSolver(
            config.linear_solver_config)

    return fluid_solver
