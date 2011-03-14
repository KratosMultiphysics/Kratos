#importing the Kratos Library
from Kratos import *
from KratosStructuralApplication import *
from KratosExternalSolversApplication import *
#from KratosMKLSolversApplication import *

def AddVariables(model_part):
    model_part.AddNodalSolutionStepVariable(DISPLACEMENT);
    model_part.AddNodalSolutionStepVariable(REACTION);
    model_part.AddNodalSolutionStepVariable(BODY_FORCE);
    model_part.AddNodalSolutionStepVariable(NEGATIVE_FACE_PRESSURE);
    model_part.AddNodalSolutionStepVariable(POSITIVE_FACE_PRESSURE);
    model_part.AddNodalSolutionStepVariable(INSITU_STRESS);
    model_part.AddNodalSolutionStepVariable(FACE_LOAD);
    model_part.AddNodalSolutionStepVariable(FORCE);
    model_part.AddNodalSolutionStepVariable(DENSITY);
    model_part.AddNodalSolutionStepVariable(ROTATION);
    model_part.AddNodalSolutionStepVariable(MOMENTUM);



    print "*********************************************************************** "
    print "Variables for the Static Structural Arc Length Solution added correctly"
    print "*********************************************************************** "    


def ChangeCondition(model_part, lamda):
	for node in model_part.Nodes:
	    new_load = node.GetSolutionStepValue(FORCE)*lamda;
	    node.SetSolutionStepValue(FORCE,0,new_load)


def AddDofs(model_part):
  for node in model_part.Nodes:
    #adding dofs
      node.AddDof(DISPLACEMENT_X,REACTION_X);
      node.AddDof(DISPLACEMENT_Y,REACTION_Y);
      node.AddDof(DISPLACEMENT_Z,REACTION_Z);
      #node.AddDof(ROTATION_X,MOMENTUM_X);
      #node.AddDof(ROTATION_Y,MOMENTUM_Y);
      #node.AddDof(ROTATION_Z,MOMENTUM_Z);
 
  print "*********************************************************************** "
  print "Dofs for the Static Structural Arc Length Solution added correctly"
  print "*********************************************************************** "

class StaticStructuralSolver:
    #######################################################################
    def __init__(self,model_part,domain_size):

        self.model_part    = model_part
        self.time_scheme   = ResidualBasedIncrementalUpdateStaticScheme()

        
        #self.time_scheme   = ParallelResidualBasedIncrementalUpdateStaticScheme()

	## Varibles de Control de Arc Lenght Method
	self.Ide                        = 5
	self.factor_delta_lmax          = 1.00
	self.toler                      = 1.0E-9
        self.norm                       = 1.0E-5 
	self.max_iter                   = 20
       
	
         #self.time_scheme = ResidualBasedPredictorCorrectorBossakScheme(self.damp_factor)
        #definition of the solvers. Super_Lu Default
        self.structure_linear_solver     =   SkylineLUFactorizationSolver()
        #self.structure_linear_solver      =   SuperLUSolver()
        #self.structure_linear_solver =  MKLPardisoSolver()
        #pDiagPrecond = ParallelDiagonalPreconditioner()
        #self.structure_linear_solver     =  ParallelCGSolver(1e-8, 5000,pDiagPrecond)
	#self.structure_linear_solver    =   Preconditioner()
        #self.structure_linear_solver    =   IterativeSolver() 

        #pDiagPrecond = DiagonalPreconditioner()
        #LST  = 1E-9
        #LSMI = 5000  
        #self.structure_linear_solver  =  BICGSTABSolver(LST,LSMI,pDiagPrecond)

        #definition of the convergence criteria
        Displacement   =   DisplacementCriteria(1E-6,1E-9)
        Residual       =   ResidualCriteria(1E-3,1E-6)

        self.conv_criteria      = AndCriteria(Residual, Displacement)
        #self.conv_Residual     = ResidualCriteria(0.000001,1E-9)
	#self.conv_Displacement = DisplacementCriteria(self.norm,self.toler)
        #self.conv_criteria     = ResDisCriteria(self.conv_Residual, self.conv_Displacement)
	#self.conv_criteria     = ResidualDisplacementCriteria(self.norm,self.toler)
	#self.conv_criteria     = DisplacementCriteria(self.norm,self.toler)
	#self.conv_criteria    = ParallelDisplacementCriteria(0.000001,1e-9)
        #self.conv_criteria     = DisplacementCriteria(self.norm,self.toler)

        #definition of the convergence criteria
       

        self.CalculateReactionFlag  = True
        self.ReformDofSetAtEachStep = True
        self.MoveMeshFlag           = True
	self.ApplyBodyForce         = False
        
    #######################################################################
    def Initialize(self):
        
        self.solver = 	ResidualBasedArcLenghtStrategy(self.model_part,self.time_scheme,self.structure_linear_solver,self.conv_criteria,self.Ide,self.max_iter,self.factor_delta_lmax, self.CalculateReactionFlag, self.ReformDofSetAtEachStep, self.MoveMeshFlag,self.ApplyBodyForce)
       
         
                 
    #######################################################################   
    def Solve(self):
        (self.solver).Solve()

    #######################################################################   
    def SetEchoLevel(self,level):
        (self.solver).SetEchoLevel(level)

    def ChangeCondition(self, model_part, lamda):
              for cond in model_part.Conditions:
                   print cond


      

