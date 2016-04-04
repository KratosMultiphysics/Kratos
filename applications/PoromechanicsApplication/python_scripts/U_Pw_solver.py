from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
#importing the Kratos Library
from KratosMultiphysics import *
from KratosMultiphysics.PoromechanicsApplication import *
from KratosMultiphysics.ExternalSolversApplication import *
#check that KratosMultiphysics was imported in the main script
CheckForPreviousImport()


def AddVariables(model_part):
    ##Solid Variables
    #add displacements
    model_part.AddNodalSolutionStepVariable(DISPLACEMENT)
    #add reactions for the displacements
    model_part.AddNodalSolutionStepVariable(REACTION)
    #add dynamic variables
    model_part.AddNodalSolutionStepVariable(VELOCITY)
    model_part.AddNodalSolutionStepVariable(ACCELERATION)
    #add variables for the solid conditions
    model_part.AddNodalSolutionStepVariable(FORCE)
    model_part.AddNodalSolutionStepVariable(FACE_LOAD)
    model_part.AddNodalSolutionStepVariable(NORMAL_CONTACT_STRESS)
    model_part.AddNodalSolutionStepVariable(TANGENTIAL_CONTACT_STRESS)
    ##Fluid Variables
    #add water pressure
    model_part.AddNodalSolutionStepVariable(WATER_PRESSURE)
    #add reactions for the water pressure
    model_part.AddNodalSolutionStepVariable(REACTION_WATER_PRESSURE)
    #add dynamic variables
    model_part.AddNodalSolutionStepVariable(DT_WATER_PRESSURE)
    #add variables for the water conditions
    model_part.AddNodalSolutionStepVariable(NORMAL_FLUID_FLUX)
    ##Common variables
    model_part.AddNodalSolutionStepVariable(VOLUME_ACCELERATION)
    
    print("Variables correctly added")


def AddDofs(model_part):
    for node in model_part.Nodes:
        ##Solid dofs
        node.AddDof(DISPLACEMENT_X,REACTION_X)
        node.AddDof(DISPLACEMENT_Y,REACTION_Y)
        node.AddDof(DISPLACEMENT_Z,REACTION_Z)
        ##Fluid dofs
        node.AddDof(WATER_PRESSURE,REACTION_WATER_PRESSURE);

    print("DOFs correctly added")


def CreateSolver(model_part, config):
    
    main_solver = UPwSolver()
    
    #Setting the linear solver (direct, iterative...)
    if(config.linear_solver == "Direct"):
        if(config.direct_solver == "Super_LU"):
            linear_solver = SuperLUSolver()
        elif(config.direct_solver == "Skyline_LU_factorization"):
            linear_solver = SkylineLUFactorizationSolver()
    elif(config.linear_solver == "Iterative"):
        if(config.iterative_solver == "BICGSTAB"):
            tolerance = 1e-5
            max_iterations = 1000
            precond = ILU0Preconditioner()
            linear_solver = BICGSTABSolver(tolerance,max_iterations,precond)
        elif(config.iterative_solver == "AMGCL"):
            tolerance = 1e-5
            max_iterations = 1000
            verbosity = 0 #0->shows no information, 1->some information, 2->all the information
            gmres_size = 50
            linear_solver = AMGCLSolver(AMGCLSmoother.ILU0,AMGCLIterativeSolverType.BICGSTAB,tolerance,max_iterations,verbosity,gmres_size)
    
    #Setting the solution scheme (quasi-static, dynamic)
    beta = 0.25
    gamma = 0.5
    theta = 0.5 #(0: Forward Euler, 1: Backward Euler, 0.5: Crank-Nicolson)
    if(config.analysis_type == "Quasi-Static"):
        solution_scheme = NewmarkQuasistaticUPwScheme(model_part,beta,gamma,theta)
    elif(config.analysis_type == "Dynamic"):
        rayleigh_m = 0.0
        rayleigh_k = 0.0
        solution_scheme = NewmarkDynamicUPwScheme(model_part,beta,gamma,theta,rayleigh_m,rayleigh_k)
    
    #Setting the builder_and_solver
    builder_and_solver = ResidualBasedEliminationBuilderAndSolver(linear_solver)
    if(config.linear_solver == "Iterative" and config.iterative_solver == "AMGCL"):
        builder_and_solver = ResidualBasedBlockBuilderAndSolver(linear_solver)
    
    #Setting the strategy type
    dofs_rel_tol = config.dofs_relative_tolerance
    residual_rel_tol = config.residual_relative_tolerance
    max_iters = config.max_iteration
    compute_react = config.compute_reactions
    if(config.strategy_type == "Newton-Raphson"):
        main_solver.strategy = NewtonRaphsonStrategy(model_part,solution_scheme,builder_and_solver,dofs_rel_tol,residual_rel_tol,max_iters,compute_react,
                                                    main_solver.reform_step_dofs,main_solver.move_mesh_flag)
    elif(config.strategy_type == "Arc-Length"):
        #Arc-length parameters
        desired_iters = 4 #The larger this number, the larger the radius
        max_rad = 20 #Times the initial radius. In order to use a constant radius, choose the same value for Min Radius and Max Radius
        min_rad = 0.5 #Times the initial radius. In order to use a constant radius, choose the same value for Min Radius and Max Radius
        main_solver.strategy = RammArcLengthStrategy(model_part,solution_scheme,builder_and_solver,dofs_rel_tol,residual_rel_tol,max_iters,desired_iters,
                                                    max_rad,min_rad,compute_react,main_solver.reform_step_dofs,main_solver.move_mesh_flag)
    
    return main_solver


class UPwSolver:

    def __init__(self):
                
        #Default level of echo for the solving strategy
         # 0 -> mute... no echo at all
         # 1 -> printing time and basic informations
         # 2 -> printing linear solver data
         # 3 -> Print of debug informations: Echo of stiffness matrix, Dx, b...
        self.echo_level = 0

        #Default strtategy options
        self.reform_step_dofs = False
        self.move_mesh_flag = True

    def Initialize(self):

        #Set echo level
        self.strategy.SetEchoLevel(self.echo_level)

        #Check if everything is assigned correctly
        self.strategy.Check()
        
        #Initialize strategy
        self.strategy.Initialize()
        
    def Solve(self):
        self.strategy.Solve()

    def Finalize(self):
        self.strategy.Clear()
