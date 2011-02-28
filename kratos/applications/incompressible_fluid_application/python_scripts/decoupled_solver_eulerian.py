#importing the Kratos Library
from Kratos import *
from KratosIncompressibleFluidApplication import *
##from KratosExternalSolversApplication import * # For SuperLU solver


def AddVariables(model_part):
    model_part.AddNodalSolutionStepVariable(VELOCITY);
    model_part.AddNodalSolutionStepVariable(ACCELERATION);
    model_part.AddNodalSolutionStepVariable(MESH_VELOCITY);
    model_part.AddNodalSolutionStepVariable(PRESSURE);
##    model_part.AddNodalSolutionStepVariable(AIR_PRESSURE);
##    model_part.AddNodalSolutionStepVariable(IS_FLUID);
##    model_part.AddNodalSolutionStepVariable(IS_POROUS);
##    model_part.AddNodalSolutionStepVariable(IS_STRUCTURE);
##    model_part.AddNodalSolutionStepVariable(IS_FREE_SURFACE);
    model_part.AddNodalSolutionStepVariable(IS_INTERFACE);
##    model_part.AddNodalSolutionStepVariable(IS_BOUNDARY);
    model_part.AddNodalSolutionStepVariable(DISPLACEMENT);
    model_part.AddNodalSolutionStepVariable(VISCOSITY);
    model_part.AddNodalSolutionStepVariable(DENSITY);
##    model_part.AddNodalSolutionStepVariable(POROSITY);
##    model_part.AddNodalSolutionStepVariable(DENSITY_AIR);
##    model_part.AddNodalSolutionStepVariable(AIR_SOUND_VELOCITY);
##    model_part.AddNodalSolutionStepVariable(SOUND_VELOCITY);
    model_part.AddNodalSolutionStepVariable(BODY_FORCE);
    model_part.AddNodalSolutionStepVariable(NODAL_AREA);
##    model_part.AddNodalSolutionStepVariable(NODAL_H);
    model_part.AddNodalSolutionStepVariable(ADVPROJ);
    model_part.AddNodalSolutionStepVariable(DIVPROJ);
##    model_part.AddNodalSolutionStepVariable(THAWONE);
##    model_part.AddNodalSolutionStepVariable(THAWTWO); 
    model_part.AddNodalSolutionStepVariable(REACTION);
    model_part.AddNodalSolutionStepVariable(REACTION_WATER_PRESSURE);
##    model_part.AddNodalSolutionStepVariable(REACTION_AIR_PRESSURE);
##    model_part.AddNodalSolutionStepVariable(EXTERNAL_PRESSURE);
##    model_part.AddNodalSolutionStepVariable(WATER_PRESSURE);
##    model_part.AddNodalSolutionStepVariable(AIR_PRESSURE_DT);
##    model_part.AddNodalSolutionStepVariable(ARRHENIUS); 
    model_part.AddNodalSolutionStepVariable(FLAG_VARIABLE);
##    model_part.AddNodalSolutionStepVariable(NORMAL);


    print "variables for the dynamic structural solution added correctly"
        
def AddDofs(model_part):
    for node in model_part.Nodes:
        #adding dofs
        node.AddDof(VELOCITY_X,REACTION_X);
        node.AddDof(VELOCITY_Y,REACTION_Y);
        node.AddDof(VELOCITY_Z,REACTION_Z);
        node.AddDof(PRESSURE,REACTION_WATER_PRESSURE);
##        node.AddDof(AIR_PRESSURE,REACTION_AIR_PRESSURE);
        
    print "dofs for the monolithic solver added correctly"

class DecoupledSolver:
    #######################################################################
    def __init__(self,model_part,domain_size):

        self.model_part = model_part

        self.alpha = -0.3
        self.move_mesh_strategy = 0
        self.time_scheme = ResidualBasedPredictorCorrectorVelocityBossakScheme\
                           ( self.alpha,self.move_mesh_strategy )
        #definition of the solvers
##        self.linear_solver =  SkylineLUFactorizationSolver()
        

        self.linear_tol = 1e-3

        vel_Precond = DiagonalPreconditioner()
##        vel_Precond = ILU0Preconditioner()
        self.velocity_linear_solver = BICGSTABSolver(self.linear_tol,\
                                                     5000,vel_Precond)

        pr_Precond = DiagonalPreconditioner()
##        pr_Precond = ILU0Preconditioner()
        self.pressure_linear_solver = BICGSTABSolver(self.linear_tol,\
                                                     5000,pr_Precond)
##        self.pressure_linear_solver =SuperLUSolver()
        
        #definition of the convergence criteria
        self.rel_vel_tol = 1e-5
        self.abs_vel_tol = 1e-7
        self.rel_pres_tol = 1e-5
        self.abs_pres_tol = 1e-7

        self.max_iter = 10

        self.dynamic_tau = 0.0
        self.oss_switch  = 0
                            
        #default settings
        self.echo_level = 0
        self.CalculateReactionFlag = True
        self.ReformDofSetAtEachStep = False
##        self.CalculateNormDxFlag = True
        self.MoveMeshFlag = False
        
        self.velocity_correction = 2
        # 0: divergence-free condition is not explicitly imposed
        # 1: divergence-free condition imposed by a diagonal divergence matrix
        # 2: divergence-free condition imposed by the full divergence operator

        # inexact Newton iterations (to use, call self.UseInexactNewtonScheme())
        self.use_inexact_newton=False
        self.IN_min_tol = self.linear_tol
        self.IN_max_tol=0.1
        self.IN_gamma=0.9
        
    #######################################################################
    def Initialize(self):
        #creating the solution strategy

        self.builder_and_solver = PressureSplittingBuilderAndSolver\
                                  (self.velocity_linear_solver,\
                                   self.pressure_linear_solver,\
                                   self.velocity_correction,\
                                   self.use_inexact_newton,\
                                   self.IN_min_tol,\
                                   self.IN_max_tol,\
                                   self.IN_gamma)

        self.conv_criteria = VelPrCriteria(self.rel_vel_tol,self.abs_vel_tol,\
                                        self.rel_pres_tol,self.abs_pres_tol)
##        self.conv_criteria = UPCriteria(self.rel_vel_tol,self.abs_vel_tol,
##                                        self.rel_pres_tol,self.abs_pres_tol)


        # Note that the strategy asks for a solver but doesn't use it (when
        # called using this constructor). This is good, as this builder and
        # solver uses 2 different solvers and one of them is given here
        # arbitrarily
        self.solver = ResidualBasedNewtonRaphsonStrategy\
                      (self.model_part,\
                       self.time_scheme,\
                       self.pressure_linear_solver,\
                       self.conv_criteria,\
                       self.builder_and_solver,\
                       self.max_iter,\
                       self.CalculateReactionFlag,\
                       self.ReformDofSetAtEachStep,\
                       self.MoveMeshFlag)
        
        (self.solver).SetEchoLevel(self.echo_level)

        self.model_part.ProcessInfo.SetValue(DYNAMIC_TAU, self.dynamic_tau);
        self.model_part.ProcessInfo.SetValue(OSS_SWITCH, self.oss_switch );
	                     
    #######################################################################   
    def Solve(self):
        self.solver.Solve()
##        self.builder_and_solver.SetUpDofSet(self.time_scheme,self.model_part)
       

    #######################################################################   
    def SetEchoLevel(self,level):
        (self.solver).SetEchoLevel(level)
    
    ########################################################################

    def UseInexactNewtonScheme(self,min_tol=-1,\
                               max_tol=0.1,gamma=0.9):
        self.use_inexact_newton=True
        if (min_tol == -1):
            self.IN_min_tol = self.linear_tol
        else:
            self.IN_min_tol=min_tol
        self.IN_max_tol=max_tol
        self.IN_gamma=gamma

    def SetRebuildLevel(self,level=1):
        self.rebuild_level = level
        if (level != 0) and (self.ReformDofSetAtEachStep == True):
            self.ReformDofSetAtEachStep = False
            print "WARNING: In DecoupledSolver.SetRebuildLevel():\n"+\
                  "I am setting ReformDofSetAtEachStep = False to be able to"+\
                  "reuse the system matrix structure. Use this option only if"+\
                  "you don't intend to reform the Dof set (for example due to"+\
                  "remeshing)"

    def UseVelocityCorrection(self,level=2):
        self.velocity_correction = level
        




