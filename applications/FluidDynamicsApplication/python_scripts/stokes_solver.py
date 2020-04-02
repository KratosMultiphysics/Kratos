# importing the Kratos Library
import KratosMultiphysics as kratoscore
import KratosMultiphysics.FluidDynamicsApplication as cfd
import KratosMultiphysics.python_linear_solver_factory as linear_solver_factory


def AddVariables(model_part, settings=None):
    model_part.AddNodalSolutionStepVariable(kratoscore.VELOCITY)
    model_part.AddNodalSolutionStepVariable(kratoscore.PRESSURE)
    #model_part.AddNodalSolutionStepVariable(kratoscore.VISCOSITY)
    model_part.AddNodalSolutionStepVariable(kratoscore.DENSITY)
    model_part.AddNodalSolutionStepVariable(kratoscore.BODY_FORCE) #TODO: decide if it is needed. if constant it could be passed in properties
    model_part.AddNodalSolutionStepVariable(kratoscore.REACTION) #in case this variable could be removed if no reactions must be computed
    #model_part.AddNodalSolutionStepVariable(kratoscore.REACTION_WATER_PRESSURE) #in case this variable could be removed if no reactions must be computed
    model_part.AddNodalSolutionStepVariable(kratoscore.EXTERNAL_PRESSURE)
    model_part.AddNodalSolutionStepVariable(kratoscore.NORMAL) #TODO: this variable is not strictly needed by the solver - may be needed by other utilities
    #model_part.AddNodalSolutionStepVariable(kratoscore.IS_STRUCTURE) #TODO: remove as deprecated!!
    #model_part.AddNodalSolutionStepVariable(kratoscore.MESH_VELOCITY) #TODO: remove. needed because of the Condition used
    #model_part.AddNodalSolutionStepVariable(kratoscore.ACCELERATION) #TODO: remove! needed because of the Condition used
    print("variables for the  monolithic solver symbolic added correctly")


def AddDofs(model_part, settings=None):
    for node in model_part.Nodes:
        # adding dofs
        node.AddDof(kratoscore.VELOCITY_X, kratoscore.REACTION_X)
        node.AddDof(kratoscore.VELOCITY_Y, kratoscore.REACTION_Y)
        node.AddDof(kratoscore.VELOCITY_Z, kratoscore.REACTION_Z)
        node.AddDof(kratoscore.PRESSURE, kratoscore.REACTION_WATER_PRESSURE)


    print("dofs for the vms monolithic solver added correctly")

def CreateSolver(model_part, settings):
    fluid_solver = StokesSolver(model_part, settings)
    return fluid_solver

class StokesSolver:

    def __init__(self, model_part, settings):

        self.model_part = model_part

        #note that all settingsuration parameters MUST be passed.
        self.domain_size = settings.domain_size
        self.rel_vel_tol = settings.vel_tolerance
        self.abs_vel_tol = settings.abs_vel_tolerance
        self.rel_pres_tol = settings.press_tolerance
        self.abs_pres_tol = settings.abs_press_tolerance
        self.dynamic_tau = settings.dynamic_tau
        self.max_iter = settings.max_iteration
        self.echo_level = settings.echo_level
        self.compute_reactions = settings.compute_reactions
        self.reform_dofs_at_each_step = settings.reform_dofs_at_each_step

        self.linear_solver = linear_solver_factory.ConstructSolver(settings.linear_solver_settings)

        time_order = 2
        self.time_discretization = KratosMultiphysics.TimeDiscretization.BDF(time_order)

        self.conv_criteria = cfd.VelPrCriteria(self.rel_vel_tol, self.abs_vel_tol,
                                               self.rel_pres_tol, self.abs_pres_tol)

        (self.conv_criteria).SetEchoLevel(self.echo_level)

        self.time_scheme = kratoscore.ResidualBasedIncrementalUpdateStaticScheme()

        builder_and_solver = kratoscore.ResidualBasedBlockBuilderAndSolver(self.linear_solver)

        move_mesh_flag = False #user should NOT configure this
        self.fluid_solver = kratoscore.ResidualBasedNewtonRaphsonStrategy(
            self.model_part, self.time_scheme, self.linear_solver, self.conv_criteria,
            builder_and_solver, self.max_iter, self.compute_reactions, self.reform_dofs_at_each_step, move_mesh_flag)
        (self.fluid_solver).SetEchoLevel(self.echo_level)
        self.fluid_solver.Check()

        self.model_part.ProcessInfo.SetValue(kratoscore.DYNAMIC_TAU, self.dynamic_tau)

        print("Construction stokes solver finished")

    #
    def Initialize(self):
        print ("Initialization stokes solver finished")

    def Solve(self):
        (self.time_discretization).ComputeAndSaveBDFCoefficients(self.model_part.ProcessInfo)
        self.fluid_solver.Solve()

    def SetEchoLevel(self, level):
        self.fluid_solver.SetEchoLevel(level)

    def Clear(self):
        self.fluid_solver.Clear()

    def Check(self):
        self.fluid_solver.Check()
