from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# importing the Kratos Library
from KratosMultiphysics import *
from KratosMultiphysics.IncompressibleFluidApplication import *
from KratosMultiphysics.FluidDynamicsApplication import *
# Check that KratosMultiphysics was imported in the main script
CheckForPreviousImport()


def AddVariables(model_part, settings=None):
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


def AddDofs(model_part, settings=None):
    for node in model_part.Nodes:
        # adding dofs
        node.AddDof(VELOCITY_X, REACTION_X)
        node.AddDof(VELOCITY_Y, REACTION_Y)
        node.AddDof(VELOCITY_Z, REACTION_Z)
        node.AddDof(PRESSURE, REACTION_WATER_PRESSURE)


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
        
        import linear_solver_factory
        self.linear_solver = linear_solver_factory.ConstructSolver(settings.linear_solver_settings)
        
        self.bdf_process = ComputeBDFCoefficientsProcess(model_part,2)
        
        self.conv_criteria = VelPrCriteria(self.rel_vel_tol, self.abs_vel_tol,
                                           self.rel_pres_tol, self.abs_pres_tol)
            
        self.time_scheme = ResidualBasedIncrementalUpdateStaticScheme() 
        
        builder_and_solver = ResidualBasedBlockBuilderAndSolver(self.linear_solver)
        
        move_mesh_flag = False #user should NOT configure this
        self.fluid_solver = ResidualBasedNewtonRaphsonStrategy(
            self.model_part, self.time_scheme, self.linear_solver, self.conv_criteria,
            builder_and_solver, self.max_iter, self.compute_reactions, self.reform_dofs_at_each_step, move_mesh_flag)
        (self.fluid_solver).SetEchoLevel(self.echo_level)
        self.fluid_solver.Check()

        self.model_part.ProcessInfo.SetValue(DYNAMIC_TAU, self.dynamic_tau)

        print("Construction stokes solver finished")

    #
    def Initialize(self):
        print ("Initialization stokes solver finished")
    
    def Solve(self):
        self.bdf_process.Execute()
        self.fluid_solver.Solve()

    def SetEchoLevel(self, level):
        self.fluid_solver.SetEchoLevel(level)

    def Clear(self):
        self.fluid_solver.Clear()
        
    def Check(self):
        self.fluid_solver.Check()


