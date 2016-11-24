from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
#importing the Kratos Library
from KratosMultiphysics import *
from KratosMultiphysics.ExternalSolversApplication import *
from KratosMultiphysics.SolidMechanicsApplication import *
from KratosMultiphysics.DamApplication import *
from KratosMultiphysics.PoromechanicsApplication import *
#check that KratosMultiphysics was imported in the main script
CheckForPreviousImport()


def AddVariables(model_part):
    #add displacements
    model_part.AddNodalSolutionStepVariable(DISPLACEMENT)
    #add reactions for the displacements
    model_part.AddNodalSolutionStepVariable(REACTION)
    #add dynamic variables
    model_part.AddNodalSolutionStepVariable(VELOCITY)
    model_part.AddNodalSolutionStepVariable(ACCELERATION)
    #add variables for the solid conditions
    model_part.AddNodalSolutionStepVariable(POINT_LOAD)
    model_part.AddNodalSolutionStepVariable(LINE_LOAD)
    model_part.AddNodalSolutionStepVariable(SURFACE_LOAD)
    model_part.AddNodalSolutionStepVariable(NEGATIVE_FACE_PRESSURE)
    #add volume acceleration
    model_part.AddNodalSolutionStepVariable(VOLUME_ACCELERATION)
    
    model_part.AddNodalSolutionStepVariable(NODAL_AREA)
    model_part.AddNodalSolutionStepVariable(NODAL_CAUCHY_STRESS_TENSOR)
    model_part.AddNodalSolutionStepVariable(Vi_POSITIVE)
    model_part.AddNodalSolutionStepVariable(Viii_POSITIVE)
    model_part.AddNodalSolutionStepVariable(NODAL_JOINT_WIDTH)
    model_part.AddNodalSolutionStepVariable(NODAL_JOINT_AREA)

    print("Variables correctly added")


def AddDofs(model_part):
    for node in model_part.Nodes:
        ##Solid dofs
        node.AddDof(DISPLACEMENT_X,REACTION_X)
        node.AddDof(DISPLACEMENT_Y,REACTION_Y)
        node.AddDof(DISPLACEMENT_Z,REACTION_Z)

    print("DOFs correctly added")


def CreateSolver(model_part, config, results):

    dam_mechanical_solver = DamMechanicalSolver()
    
    #Setting the strategy type
    convergence_criterion = config.convergence_criterion
    dis_rel_tol = config.displacement_rel_tol
    dis_abs_tol = config.displacement_abs_tol
    res_rel_tol = config.residual_rel_tol
    res_abs_tol = config.residual_abs_tol
    echo_level = config.echo_level
        
    max_iters = config.max_iteration
    compute_react = config.compute_reactions

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

    #Setting the solution scheme (quasi-static, dynamic) and the possibility to plot Nodal Cauchy stress tensor
    if("NODAL_CAUCHY_STRESS_TENSOR" in results):
        if (config.analysis_type == "Quasi-Static"):       
            solution_scheme = IncrementalUpdateStaticSmoothingScheme()
        else:
            damp_factor_m = -0.01
            solution_scheme = BossakDisplacementSmoothingScheme(damp_factor_m)
    else:
        if (config.analysis_type == "Quasi-Static"):
            solution_scheme = ResidualBasedIncrementalUpdateStaticScheme()
        else:
            damp_factor_m = -0.01
            solution_scheme = ResidualBasedBossakDisplacementScheme(damp_factor_m)     
        
    #Setting the builder_and_solver
    builder_and_solver = ResidualBasedEliminationBuilderAndSolver(linear_solver)
    if(config.linear_solver == "Iterative" and config.iterative_solver == "AMGCL"):
        builder_and_solver = ResidualBasedBlockBuilderAndSolver(linear_solver)        

    #Convergence Criterion
    if(convergence_criterion == "Displacement_criterion"):
        convergence_criterion = DisplacementConvergenceCriterion(dis_rel_tol, dis_abs_tol)
        convergence_criterion.SetEchoLevel(echo_level)
    elif(convergence_criterion == "Residual_criterion"):
        convergence_criterion = ResidualCriteria(res_rel_tol, res_abs_tol)
        convergence_criterion.SetEchoLevel(echo_level)
    elif(convergence_criterion == "And_criterion"):
        Displacement = DisplacementConvergenceCriterion(dis_rel_tol, dis_abs_tol)
        Displacement.SetEchoLevel(echo_level)
        Residual = ResidualCriteria(res_rel_tol, res_abs_tol)
        Residual.SetEchoLevel(echo_level)
        convergence_criterion = AndCriteria(Residual, Displacement)
        
    if(config.strategy_type == "Newton-Raphson"):
        
        model_part.ProcessInfo[IS_CONVERGED]=True
        dam_mechanical_solver.strategy = ResidualBasedNewtonRaphsonStrategy(model_part,solution_scheme,linear_solver,convergence_criterion,builder_and_solver,max_iters,
                                                                            compute_react,dam_mechanical_solver.reform_step_dofs,dam_mechanical_solver.move_mesh_flag)


    return dam_mechanical_solver


class DamMechanicalSolver:

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
