from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# importing the Kratos Library
from KratosMultiphysics import *
from KratosMultiphysics.IncompressibleFluidApplication import *
from KratosMultiphysics.FluidDynamicsApplication import *
# Check that KratosMultiphysics was imported in the main script
CheckForPreviousImport()


def AddVariables(model_part):
    model_part.AddNodalSolutionStepVariable(VELOCITY)
    model_part.AddNodalSolutionStepVariable(ACCELERATION)
    model_part.AddNodalSolutionStepVariable(MESH_VELOCITY)
    model_part.AddNodalSolutionStepVariable(PRESSURE)
    model_part.AddNodalSolutionStepVariable(AIR_PRESSURE)
    model_part.AddNodalSolutionStepVariable(IS_FLUID)
    model_part.AddNodalSolutionStepVariable(IS_POROUS)
    model_part.AddNodalSolutionStepVariable(IS_STRUCTURE)
    model_part.AddNodalSolutionStepVariable(IS_FREE_SURFACE)
    model_part.AddNodalSolutionStepVariable(IS_INTERFACE)
    model_part.AddNodalSolutionStepVariable(IS_BOUNDARY)
    model_part.AddNodalSolutionStepVariable(DISPLACEMENT)
    model_part.AddNodalSolutionStepVariable(VISCOSITY)
    model_part.AddNodalSolutionStepVariable(DENSITY)
    model_part.AddNodalSolutionStepVariable(POROSITY)
    model_part.AddNodalSolutionStepVariable(DENSITY_AIR)
    model_part.AddNodalSolutionStepVariable(AIR_SOUND_VELOCITY)
    model_part.AddNodalSolutionStepVariable(SOUND_VELOCITY)
    model_part.AddNodalSolutionStepVariable(BODY_FORCE)
    model_part.AddNodalSolutionStepVariable(NODAL_AREA)
    model_part.AddNodalSolutionStepVariable(NODAL_H)
    model_part.AddNodalSolutionStepVariable(ADVPROJ)
    model_part.AddNodalSolutionStepVariable(DIVPROJ)
    model_part.AddNodalSolutionStepVariable(THAWONE)
    model_part.AddNodalSolutionStepVariable(THAWTWO)
    model_part.AddNodalSolutionStepVariable(REACTION)
    model_part.AddNodalSolutionStepVariable(REACTION_WATER_PRESSURE)
    model_part.AddNodalSolutionStepVariable(EXTERNAL_PRESSURE)
    model_part.AddNodalSolutionStepVariable(WATER_PRESSURE)
    model_part.AddNodalSolutionStepVariable(AIR_PRESSURE_DT)
    model_part.AddNodalSolutionStepVariable(ARRHENIUS)
    model_part.AddNodalSolutionStepVariable(FLAG_VARIABLE)
    model_part.AddNodalSolutionStepVariable(NORMAL)

    print("variables for the MONOLITHIC_SOLVER_EULERIAN added correctly")


def AddDofs(model_part):
    for node in model_part.Nodes:
        # adding dofs
        node.AddDof(VELOCITY_X, REACTION_X)
        node.AddDof(VELOCITY_Y, REACTION_Y)
        node.AddDof(VELOCITY_Z, REACTION_Z)
        node.AddDof(PRESSURE, REACTION_WATER_PRESSURE)
        node.AddDof(AIR_PRESSURE, REACTION_AIR_PRESSURE)

    print("dofs for the monolithic solver added correctly")


class MonolithicSolver:
    #

    def __init__(self, model_part, domain_size):

        self.model_part = model_part

        self.alpha = -0.3
        self.move_mesh_strategy = 0
# self.time_scheme = ResidualBasedPredictorCorrectorVelocityBossakScheme( self.alpha,self.move_mesh_strategy )
        # definition of the solvers
# self.linear_solver =  SkylineLUFactorizationSolver()
# self.linear_solver =SuperLUSolver()

        pPrecond = DiagonalPreconditioner()
# pPrecond = ILU0Preconditioner()
        self.linear_solver = BICGSTABSolver(1e-6, 5000, pPrecond)

        # definition of the convergence criteria
        self.rel_vel_tol = 1e-5
        self.abs_vel_tol = 1e-7
        self.rel_pres_tol = 1e-5
        self.abs_pres_tol = 1e-7

       # self.conv_criteria = UPCriteria(1e-12,1e-14,1e-15,1e-17)
        # self.model_part.ProcessInfo.SetValue(DYNAMIC_TAU, 0.001);

        self.dynamic_tau = 0.0
        self.oss_switch = 0

        # non newtonian setting
        self.regularization_coef = 1000

        self.max_iter = 30

        # default settings
        self.echo_level = 0
        self.CalculateReactionFlag = True
        self.ReformDofSetAtEachStep = False
        self.CalculateNormDxFlag = True
        self.MoveMeshFlag = False

        # For Spalart-Allmaras
        self.turbulence_model = None
        self.domain_size = domain_size

# print "Construction monolithic solver finished"

    #
    def Initialize(self):

        # time scheme
        if self.turbulence_model is None:
            self.time_scheme = ResidualBasedPredictorCorrectorVelocityBossakSchemeTurbulent\
                (self.alpha, self.move_mesh_strategy, self.domain_size)
        else:
            self.time_scheme = ResidualBasedPredictorCorrectorVelocityBossakSchemeTurbulent\
                (self.alpha,
                 self.move_mesh_strategy,
                 self.domain_size,
                 self.turbulence_model)

        # creating the solution strategy
        self.conv_criteria = VelPrCriteria(self.rel_vel_tol, self.abs_vel_tol,
                                           self.rel_pres_tol, self.abs_pres_tol)
# self.conv_criteria = UPCriteria(self.rel_vel_tol,self.abs_vel_tol,
# self.rel_pres_tol,self.abs_pres_tol)
        self.solver = ResidualBasedNewtonRaphsonStrategy(
            self.model_part,
            self.time_scheme,
            self.linear_solver,
            self.conv_criteria,
            self.max_iter,
            self.CalculateReactionFlag,
            self.ReformDofSetAtEachStep,
            self.MoveMeshFlag)
        (self.solver).SetEchoLevel(self.echo_level)

        self.model_part.ProcessInfo.SetValue(DYNAMIC_TAU, self.dynamic_tau)
        self.model_part.ProcessInfo.SetValue(OSS_SWITCH, self.oss_switch)
        self.model_part.ProcessInfo.SetValue(M, self.regularization_coef)

# print "Initialization monolithic solver finished"

    #
    def Solve(self):
# print "*****************entering solve?????????????"
        (self.solver).Solve()
# print "solving step monolithic solver finished"

    #
    def SetEchoLevel(self, level):
        (self.solver).SetEchoLevel(level)

    #

    def ActivateSmagorinsky(self, c):
        for elem in self.model_part.Elements:
            elem.SetValue(C_SMAGORINSKY, c)

    def ActivateSpalartAllmaras(self, wall_nodes, DES=False, CDES=1.0):

        number_of_avg_elems = 10
        number_of_avg_nodes = 10
        neighbour_search = FindNodalNeighboursProcess(
            self.model_part,
            number_of_avg_elems,
            number_of_avg_nodes)
        neighbour_search.Execute()

        for node in wall_nodes:
            node.SetValue(IS_VISITED, 1.0)
            node.SetSolutionStepValue(DISTANCE, 0, 0.0)

        # Compute distance function
        distance_calculator = BodyDistanceCalculationUtils()
        distance_calculator.CalculateDistances2D(
            self.model_part.Elements, DISTANCE, 100.0)

        non_linear_tol = 0.001
        max_it = 10
        reform_dofset = self.ReformDofSetAtEachStep
        time_order = 2
        pPrecond = DiagonalPreconditioner()
        turbulence_linear_solver = BICGSTABSolver(1e-9, 5000, pPrecond)

        self.turbulence_model = SpalartAllmarasTurbulenceModel(
            self.model_part,
            turbulence_linear_solver,
            self.domain_size,
            non_linear_tol,
            max_it,
            reform_dofset,
            time_order)
        if DES:
            self.turbulence_model.ActivateDES(CDES)
