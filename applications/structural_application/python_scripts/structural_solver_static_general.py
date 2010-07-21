#importing the Kratos Library
from Kratos import *
from KratosStructuralApplication import *
from KratosExternalSolversApplication import *
from KratosMKLSolversApplication import *

def AddVariables(model_part):
    model_part.AddNodalSolutionStepVariable(DISPLACEMENT);
    model_part.AddNodalSolutionStepVariable(DISPLACEMENT_OLD);
    model_part.AddNodalSolutionStepVariable(DISPLACEMENT_NULL);
    model_part.AddNodalSolutionStepVariable(DISPLACEMENT_EINS);
    model_part.AddNodalSolutionStepVariable(DISPLACEMENT_DT);
    model_part.AddNodalSolutionStepVariable(DISPLACEMENT_NULL_DT);
    model_part.AddNodalSolutionStepVariable(DISPLACEMENT_EINS_DT);
    model_part.AddNodalSolutionStepVariable(ACCELERATION);
    model_part.AddNodalSolutionStepVariable(ACCELERATION_NULL);
    model_part.AddNodalSolutionStepVariable(ACCELERATION_EINS);
    model_part.AddNodalSolutionStepVariable(VELOCITY);
    model_part.AddNodalSolutionStepVariable(ACCELERATION);
    model_part.AddNodalSolutionStepVariable(REACTION);
    #model_part.AddNodalSolutionStepVariable(BODY_FORCE);
    model_part.AddNodalSolutionStepVariable(NEGATIVE_FACE_PRESSURE);
    model_part.AddNodalSolutionStepVariable(POSITIVE_FACE_PRESSURE);
    model_part.AddNodalSolutionStepVariable(INSITU_STRESS);
    model_part.AddNodalSolutionStepVariable(FACE_LOAD);

    print "variables for the dynamic structural solution added correctly"
        
def AddDofs(model_part):
    for node in model_part.Nodes:
        #adding dofs
        node.AddDof(DISPLACEMENT_X,REACTION_X);
        node.AddDof(DISPLACEMENT_Y,REACTION_Y);
        node.AddDof(DISPLACEMENT_Z,REACTION_Z);
    print "dofs for the dynamic structural solution added correctly"

class StaticStructuralSolver:
    #######################################################################
    def __init__(self,model_part,domain_size):

        self.model_part = model_part
        self.time_scheme = ResidualBasedIncrementalUpdateStaticScheme()

	# definition of parameters
	self.MaxLineSearchIterations    = 20
	self.MaxNewtonRapshonIterations = 500
	self.tolls     = 0.8           # energy tolerance factor on LineSearch (0.8 is ok) 
	self.amp       = 1.618         # maximum amplification factor
	self.etmxa     = 3.0           # maximum allowed step length
 	self.etmna     = 0.1           # minimum allowed step length
	self.toler     = 1.0E-9
        self.norm      = 1.0E-6

        #definition of the solvers
        #self.structure_linear_solver  =  SkylineLUFactorizationSolver()
	self.structure_linear_solver = MKLPardisoSolver() # SuperLUSolver()
        
        #definition of the convergence criteria

        #self.conv_Residual     = ResidualCriteria(0.000001,1E-9)
	#self.conv_Displacement = DisplacementCriteria(self.norm,self.toler)
        #self.conv_criteria     = ResDisCriteria(self.conv_Residual, self.conv_Displacement)
	#self.conv_criteria     = ResidualDisplacementCriteria(self.norm,self.toler)
        self.conv_criteria      = DisplacementCriteria(self.norm,self.toler)

        self.CalculateReactionFlag  = True
        self.ReformDofSetAtEachStep = True
        self.MoveMeshFlag 	    = True
        self.ApplyLineSearches      = True
	  

        
    #######################################################################
    def Initialize(self):
        #creating the solution strategy
        
        #import strategy_python
        #self.solver = strategy_python.SolvingStrategyPython(self.model_part,self.time_scheme,self.structure_linear_solver,self.conv_criteria,self.CalculateReactionFlag,self.ReformDofSetAtEachStep,self.MoveMeshFlag)
        #(self.solver).SetEchoLevel(3)

        #creating the solution strategy
        #self.solver = ResidualBasedNewtonRaphsonStrategy(self.model_part,self.time_scheme,self.structure_linear_solver,self.conv_criteria,60,True,True,True)

	self.solver = ResidualBasedNewtonRaphsonLineSearchesStrategy(self.model_part,self.time_scheme,self.structure_linear_solver,self.conv_criteria,
	self.MaxNewtonRapshonIterations,self.MaxLineSearchIterations, self.tolls, self.amp, self.etmxa, self.etmna,
        self.CalculateReactionFlag,
        self.ReformDofSetAtEachStep,
        self.MoveMeshFlag,
        self.ApplyLineSearches)        
        #(self.solver).SetReformDofSetAtEachStepFlag(True)
        #(self.solver).SetMoveMeshFlag(True)
        
 
                 
    #######################################################################   
    def Solve(self):
        (self.solver).Solve()

    #######################################################################   
    def SetEchoLevel(self,level):
        (self.solver).SetEchoLevel(level)

      

