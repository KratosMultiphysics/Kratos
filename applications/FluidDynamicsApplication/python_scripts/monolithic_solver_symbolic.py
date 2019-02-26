from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# importing the Kratos Library
from KratosMultiphysics import *
from KratosMultiphysics.IncompressibleFluidApplication import *
from KratosMultiphysics.FluidDynamicsApplication import *

def AddVariables(model_part, config=None):
    model_part.AddNodalSolutionStepVariable(VELOCITY)
    model_part.AddNodalSolutionStepVariable(PRESSURE)
    model_part.AddNodalSolutionStepVariable(VISCOSITY)
    model_part.AddNodalSolutionStepVariable(DENSITY)
    model_part.AddNodalSolutionStepVariable(BODY_FORCE)
    model_part.AddNodalSolutionStepVariable(REACTION)
    model_part.AddNodalSolutionStepVariable(REACTION_WATER_PRESSURE)
    model_part.AddNodalSolutionStepVariable(EXTERNAL_PRESSURE)
    model_part.AddNodalSolutionStepVariable(NORMAL)
    model_part.AddNodalSolutionStepVariable(Y_WALL) #TODO: remove, only needed for the wall law, and should be passed in Properties insteead
    model_part.AddNodalSolutionStepVariable(IS_STRUCTURE) #TODO: remove as deprecated!!
    model_part.AddNodalSolutionStepVariable(MESH_VELOCITY) #TODO: include in the model
    model_part.AddNodalSolutionStepVariable(ACCELERATION) #TODO: remove! needed because of the face
    print("variables for the  monolithic solver symbolic added correctly")


def AddDofs(model_part, config=None):
    for node in model_part.Nodes:
        # adding dofs
        node.AddDof(VELOCITY_X, REACTION_X)
        node.AddDof(VELOCITY_Y, REACTION_Y)
        node.AddDof(VELOCITY_Z, REACTION_Z)
        node.AddDof(PRESSURE, REACTION_WATER_PRESSURE)


    print("dofs for the vms monolithic solver added correctly")


class MonolithicSolver:

    def __init__(self, model_part, domain_size, periodic=False):

        self.model_part = model_part
        self.domain_size = domain_size

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
        self.max_iter = 5

        # default settings
        self.echo_level = 0
        self.compute_reactions = False
        self.ReformDofSetAtEachStep = False
        self.CalculateNormDxFlag = False
        self.MoveMeshFlag = False
        self.use_slip_conditions = False

        self.divergence_clearance_steps = 0

        self.bdf_process = ComputeBDFCoefficientsProcess(model_part,2)

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

        # creating the solution strategy
        self.conv_criteria = VelPrCriteria(self.rel_vel_tol, self.abs_vel_tol,
                                           self.rel_pres_tol, self.abs_pres_tol)

        self.conv_criteria.SetEchoLevel(self.echo_level)

        self.time_scheme = ResidualBasedIncrementalUpdateStaticScheme()

        builder_and_solver = ResidualBasedBlockBuilderAndSolver(self.linear_solver)

        self.fluid_solver = ResidualBasedNewtonRaphsonStrategy(
            self.model_part, self.time_scheme, self.linear_solver, self.conv_criteria,
            builder_and_solver, self.max_iter, self.compute_reactions, self.ReformDofSetAtEachStep, self.MoveMeshFlag)
        (self.fluid_solver).SetEchoLevel(self.echo_level)
        self.fluid_solver.Check()

        self.model_part.ProcessInfo.SetValue(DYNAMIC_TAU, self.dynamic_tau)


        print ("Initialization monolithic solver finished")

    def Solve(self):

        self.bdf_process.Execute()
        if self.divergence_clearance_steps > 0:
            print("Calculating divergence-free initial condition")
            # initialize with a Stokes solution step
            raise Error("yet to be implemented")

            self.divergence_clearance_steps = 0
            print("Finished divergence clearance")

        if self.ReformDofSetAtEachStep:
            if self.use_slip_conditions:
                self.normal_util.CalculateOnSimplex(
                    self.model_part, self.domain_size, IS_STRUCTURE)

        self.fluid_solver.Solve()

    #
    def SetEchoLevel(self, level):
        self.fluid_solver.SetEchoLevel(level)

    #
    def Clear(self):
        self.fluid_solver.Clear()

    def Check(self):
        self.fluid_solver.Check()



#
def CreateSolver(model_part, config, periodic=False):
    fluid_solver = MonolithicSolver(model_part, config.domain_size, periodic)

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

    import KratosMultiphysics.python_linear_solver_factory as linear_solver_factory
    if(hasattr(config, "linear_solver_config")):
        fluid_solver.linear_solver = linear_solver_factory.ConstructSolver(
            config.linear_solver_config)

    return fluid_solver
