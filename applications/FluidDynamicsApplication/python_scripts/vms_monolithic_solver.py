# importing the Kratos Library
from KratosMultiphysics import *
from KratosMultiphysics.IncompressibleFluidApplication import *
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

    print "variables for the MONOLITHIC_SOLVER_EULERIAN added correctly"


def AddDofs(model_part, config=None):
    for node in model_part.Nodes:
        # adding dofs
        node.AddDof(VELOCITY_X, REACTION_X)
        node.AddDof(VELOCITY_Y, REACTION_Y)
        node.AddDof(VELOCITY_Z, REACTION_Z)
        node.AddDof(PRESSURE, REACTION_WATER_PRESSURE)

    print "dofs for the monolithic solver added correctly"


class MonolithicSolver:
    #

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
# self.linear_solver =SuperLUSolver()
# self.linear_solver = MKLPardisoSolver()

        # pPrecond = DiagonalPreconditioner()
# pPrecond = ILU0Preconditioner()
        # self.linear_solver =  BICGSTABSolver(1e-6, 5000,pPrecond)

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

        self.turbulence_model = None

	print "Construction monolithic solver finished"

    #
    def Initialize(self):
        # check if slip conditions are defined
        if self.use_slip_conditions == False:
            for cond in self.model_part.Conditions:
                if cond.GetValue(IS_STRUCTURE) != 0.0:
                    self.use_slip_conditions = True
                    break

        # if we use slip conditions, calculate normals on the boundary
        if self.use_slip_conditions == True:
            self.normal_util = NormalCalculationUtils()
            self.normal_util.CalculateOnSimplex(
                self.model_part, self.domain_size, IS_STRUCTURE)

            for cond in self.model_part.Conditions:
                if cond.GetValue(IS_STRUCTURE) != 0.0:
                    for node in cond.GetNodes():
                        node.SetValue(IS_STRUCTURE, 1.0)

        # creating the solution strategy
        self.conv_criteria = VelPrCriteria(self.rel_vel_tol, self.abs_vel_tol,
                                           self.rel_pres_tol, self.abs_pres_tol)
# self.conv_criteria = UPCriteria(self.rel_vel_tol,self.abs_vel_tol,
# self.rel_pres_tol,self.abs_pres_tol)
        if self.turbulence_model == None:
            self.time_scheme = ResidualBasedPredictorCorrectorVelocityBossakSchemeTurbulent\
                (self.alpha, self.move_mesh_strategy,
                 self.domain_size)
        else:
            self.time_scheme = ResidualBasedPredictorCorrectorVelocityBossakSchemeTurbulent\
                (self.alpha,
                 self.move_mesh_strategy,
                 self.domain_size,
                 self.turbulence_model)

        builder_and_solver = ResidualBasedBlockBuilderAndSolver(
            self.linear_solver)

        self.solver = ResidualBasedNewtonRaphsonStrategy(
            self.model_part, self.time_scheme, self.linear_solver, self.conv_criteria,
            builder_and_solver, self.max_iter, self.compute_reactions, self.ReformDofSetAtEachStep, self.MoveMeshFlag)
        (self.solver).SetEchoLevel(self.echo_level)
        self.solver.Check()

        self.model_part.ProcessInfo.SetValue(DYNAMIC_TAU, self.dynamic_tau)
        self.model_part.ProcessInfo.SetValue(OSS_SWITCH, self.oss_switch)
        self.model_part.ProcessInfo.SetValue(M, self.regularization_coef)

# print "Initialization monolithic solver finished"
    #
    def Solve(self):
        if(self.ReformDofSetAtEachStep == True):
            if self.use_slip_conditions == True:
                self.normal_util.CalculateOnSimplex(
                    self.model_part, self.domain_size, IS_STRUCTURE)

        self.model_part.ProcessInfo.SetValue(DYNAMIC_TAU, self.dynamic_tau)
        self.model_part.ProcessInfo.SetValue(OSS_SWITCH, self.oss_switch)

        (self.solver).Solve()

    #
    def SetEchoLevel(self, level):
        (self.solver).SetEchoLevel(level)

    #
    def Clear(self):
        (self.solver).Clear()

    #
    def ActivateSmagorinsky(self, C):
        for elem in self.model_part.Elements:
            elem.SetValue(C_SMAGORINSKY, C)

    #
    def ActivateSpalartAllmaras(self, wall_nodes, DES=False, CDES=1.0):

        number_of_avg_elems = 10
        number_of_avg_nodes = 10
        neighbour_search = FindNodalNeighboursProcess(
            self.model_part, number_of_avg_elems, number_of_avg_nodes)
        neighbour_search.Execute()

        for node in wall_nodes:
            node.SetValue(IS_VISITED, 1.0)
            node.SetSolutionStepValue(DISTANCE, 0, 0.0)

        # Compute distance function
        distance_calculator = BodyDistanceCalculationUtils()
        if(self.domain_size == 2):
            distance_calculator.CalculateDistances2D(
                self.model_part.Elements, DISTANCE, 100.0)
        elif(self.domain_size == 3):
            distance_calculator.CalculateDistances3D(
                self.model_part.Elements, DISTANCE, 100.0)

        non_linear_tol = 0.001
        max_it = 10
        reform_dofset = self.ReformDofSetAtEachStep
        time_order = 2
        pPrecond = DiagonalPreconditioner()
        turbulence_linear_solver = BICGSTABSolver(1e-9, 5000, pPrecond)

        self.turbulence_model = SpalartAllmarasTurbulenceModel(
            self.model_part, turbulence_linear_solver, self.domain_size, non_linear_tol, max_it, reform_dofset, time_order)
        if DES:
            self.turbulence_model.ActivateDES(CDES)

            
       
#################################################################################################
################################################################################################# 
def CreateSolver( model_part, config ):
    fluid_solver = MonolithicSolver( model_part, config.domain_size )
    
    if( hasattr(config,"alpha") ): fluid_solver.alpha = config.alpha
    
    # definition of the convergence criteria
    if( hasattr(config,"velocity_relative_tolerance") ): fluid_solver.rel_vel_tol = config.velocity_relative_tolerance
    if( hasattr(config,"velocity_absolute_tolerance") ): fluid_solver.abs_vel_tol = config.velocity_absolute_tolerance
    if( hasattr(config,"pressure_relative_tolerance") ): fluid_solver.rel_pres_tol = config.pressure_relative_tolerance
    if( hasattr(config,"pressure_absolute_tolerance") ): fluid_solver.abs_pres_tol = config.pressure_absolute_tolerance
    if( hasattr(config,"dynamic_tau") ): fluid_solver.dynamic_tau = config.dynamic_tau
    if( hasattr(config,"oss_switch") ): fluid_solver.oss_switch = config.oss_switch
    if( hasattr(config,"max_iteration") ): fluid_solver.max_iter = config.max_iteration
    if( hasattr(config,"echo_level") ): fluid_solver.echo_level = config.echo_level
    if( hasattr(config,"compute_reactions") ): fluid_solver.compute_reactions = config.compute_reactions
    if( hasattr(config,"ReformDofSetAtEachStep") ): fluid_solver.ReformDofSetAtEachStep = config.ReformDofSetAtEachStep
    if( hasattr(config,"use_spalart_allmaras") ): fluid_solver.use_spalart_allmaras = config.use_spalart_allmaras
    if( hasattr(config,"use_des") ): fluid_solver.use_des = config.use_des
        
    import linear_solver_factory
    if( hasattr(config,"linear_solver_config") ): fluid_solver.linear_solver =  linear_solver_factory.ConstructSolver(config.linear_solver_config)
    
    return fluid_solver
