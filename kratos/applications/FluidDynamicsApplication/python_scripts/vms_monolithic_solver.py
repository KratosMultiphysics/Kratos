#importing the Kratos Library
from KratosMultiphysics import *
from KratosMultiphysics.IncompressibleFluidApplication import *
from KratosMultiphysics.FluidDynamicsApplication import *
# Check that KratosMultiphysics was imported in the main script
CheckForPreviousImport()


def AddVariables(model_part):
    model_part.AddNodalSolutionStepVariable(VELOCITY);
    model_part.AddNodalSolutionStepVariable(ACCELERATION);
    model_part.AddNodalSolutionStepVariable(MESH_VELOCITY);
    model_part.AddNodalSolutionStepVariable(PRESSURE);
    model_part.AddNodalSolutionStepVariable(IS_STRUCTURE);
    model_part.AddNodalSolutionStepVariable(DISPLACEMENT);
    model_part.AddNodalSolutionStepVariable(VISCOSITY);
    model_part.AddNodalSolutionStepVariable(DENSITY);
    model_part.AddNodalSolutionStepVariable(BODY_FORCE);
    model_part.AddNodalSolutionStepVariable(NODAL_AREA);
    model_part.AddNodalSolutionStepVariable(NODAL_H);
    model_part.AddNodalSolutionStepVariable(ADVPROJ);
    model_part.AddNodalSolutionStepVariable(DIVPROJ);
    model_part.AddNodalSolutionStepVariable(REACTION); 
    model_part.AddNodalSolutionStepVariable(REACTION_WATER_PRESSURE);
    model_part.AddNodalSolutionStepVariable(EXTERNAL_PRESSURE);
    model_part.AddNodalSolutionStepVariable(FLAG_VARIABLE);
    model_part.AddNodalSolutionStepVariable(NORMAL);
    model_part.AddNodalSolutionStepVariable(Y_WALL);


    print "variables for the MONOLITHIC_SOLVER_EULERIAN added correctly"
        
def AddDofs(model_part):
    for node in model_part.Nodes:
        #adding dofs
        node.AddDof(VELOCITY_X,REACTION_X);
        node.AddDof(VELOCITY_Y,REACTION_Y);
        node.AddDof(VELOCITY_Z,REACTION_Z);
        node.AddDof(PRESSURE,REACTION_WATER_PRESSURE);
        
    print "dofs for the monolithic solver added correctly"

class MonolithicSolver:
    #######################################################################
    def __init__(self,model_part,domain_size):

        self.model_part = model_part
        self.domain_size = domain_size

        self.alpha = -0.3
        self.move_mesh_strategy = 0
        
        #definition of the solvers
	try:
	    from KratosMultiphysics.ExternalSolversApplication import SuperLUIterativeSolver
	    self.linear_solver = SuperLUIterativeSolver()
	except:
            self.linear_solver =  SkylineLUFactorizationSolver()
##        self.linear_solver =SuperLUSolver()
##        self.linear_solver = MKLPardisoSolver()

        #pPrecond = DiagonalPreconditioner()
##        pPrecond = ILU0Preconditioner()
        #self.linear_solver =  BICGSTABSolver(1e-6, 5000,pPrecond)

        #definition of the convergence criteria
        self.rel_vel_tol = 1e-5
        self.abs_vel_tol = 1e-7
        self.rel_pres_tol = 1e-5
        self.abs_pres_tol = 1e-7

        self.dynamic_tau = 0.0
        self.oss_switch  = 0

        #non newtonian setting
        self.regularization_coef = 1000
        
        self.max_iter = 30
                            
        #default settings
        self.echo_level = 0
        self.compute_reactions = True
        self.ReformDofSetAtEachStep = True
        self.CalculateNormDxFlag = True
        self.MoveMeshFlag = False
        self.use_slip_conditions = False
    
##        print "Construction monolithic solver finished"
        
    #######################################################################
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
            self.normal_util.CalculateOnSimplex(self.model_part,self.domain_size,IS_STRUCTURE)
            
            for cond in self.model_part.Conditions:
		if cond.GetValue(IS_STRUCTURE) != 0.0:
		   for node in cond.GetNodes() :
		      node.SetValue(IS_STRUCTURE,1.0)
		   
        #creating the solution strategy
        self.conv_criteria = VelPrCriteria(self.rel_vel_tol,self.abs_vel_tol,\
                                           self.rel_pres_tol,self.abs_pres_tol)
##        self.conv_criteria = UPCriteria(self.rel_vel_tol,self.abs_vel_tol,
##                                        self.rel_pres_tol,self.abs_pres_tol)
	self.time_scheme = ResidualBasedPredictorCorrectorVelocityBossakSchemeTurbulent( self.alpha,self.move_mesh_strategy, self.domain_size )
	
        self.solver = ResidualBasedNewtonRaphsonStrategy(self.model_part,self.time_scheme,self.linear_solver,self.conv_criteria,self.max_iter,self.compute_reactions, self.ReformDofSetAtEachStep,self.MoveMeshFlag)   
        (self.solver).SetEchoLevel(self.echo_level)
        self.solver.Check()

        self.model_part.ProcessInfo.SetValue(DYNAMIC_TAU, self.dynamic_tau);
        self.model_part.ProcessInfo.SetValue(OSS_SWITCH, self.oss_switch );
        self.model_part.ProcessInfo.SetValue(M, self.regularization_coef );
        


##        print "Initialization monolithic solver finished"
	                     
    #######################################################################   
    def Solve(self):
	if(self.ReformDofSetAtEachStep == True):
            if self.use_slip_conditions == True:
                self.normal_util.CalculateOnSimplex(self.model_part,self.domain_size,IS_STRUCTURE)

        self.model_part.ProcessInfo.SetValue(DYNAMIC_TAU, self.dynamic_tau);
        self.model_part.ProcessInfo.SetValue(OSS_SWITCH, self.oss_switch );

        (self.solver).Solve()
       

    #######################################################################   
    def SetEchoLevel(self,level):
        (self.solver).SetEchoLevel(level)
    
    ########################################################################
    def Clear(self):
        (self.solver).Clear()
        
    ########################################################################
    def ActivateSmagorinsky(self,C):
        for elem in self.model_part.Elements:
            elem.SetValue(C_SMAGORINSKY,C)



