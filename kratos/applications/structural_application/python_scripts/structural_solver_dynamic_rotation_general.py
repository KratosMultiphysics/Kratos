#importing the Kratos Library
from KratosMultiphysics import *
from KratosMultiphysics.StructuralApplication import *


def AddVariables(model_part):
    model_part.AddNodalSolutionStepVariable(DISPLACEMENT);
    model_part.AddNodalSolutionStepVariable(VELOCITY);
    model_part.AddNodalSolutionStepVariable(ACCELERATION);
    model_part.AddNodalSolutionStepVariable(ROTATION);
    model_part.AddNodalSolutionStepVariable(ANGULAR_VELOCITY);
    model_part.AddNodalSolutionStepVariable(ANGULAR_ACCELERATION);
    model_part.AddNodalSolutionStepVariable(REACTION);
    model_part.AddNodalSolutionStepVariable(NEGATIVE_FACE_PRESSURE);
    model_part.AddNodalSolutionStepVariable(POSITIVE_FACE_PRESSURE);


    print "variables for the dynamic structural solution added correctly"
        
def AddDofs(model_part):
    for node in model_part.Nodes:
        #adding dofs
        node.AddDof(DISPLACEMENT_X);
        node.AddDof(DISPLACEMENT_Y);
        node.AddDof(DISPLACEMENT_Z);
        node.AddDof(ROTATION_X);
        node.AddDof(ROTATION_Y);
        node.AddDof(ROTATION_Z);
    print "dofs for the dynamic structural solution added correctly"


class DynamicStructuralSolver:
    #######################################################################
    def __init__(self,model_part,domain_size):

        self.model_part = model_part
        self.echo_level = 0
        
        self.damp_factor = -0.1
        self.toll = 1e-6
        self.absolute_tol = 1e-9
       	self.MaxLineSearchIterations    = 20
	self.MaxNewtonRapshonIterations = 500
	self.tolls     = 0.8           # energy tolerance factor on LineSearch (0.8 is ok) 
	self.amp       = 1.618         # maximum amplification factor
	self.etmxa     = 3.0           # maximum allowed step length
 	self.etmna     = 0.1           # minimum allowed step length
	self.toler     = 1.0E-9
        self.norm      = 1.0E-6


        #definition of the solvers
        self.structure_linear_solver =  SkylineLUFactorizationSolver()

        self.CalculateReactionFlag  = True
        self.ReformDofSetAtEachStep = True
        self.MoveMeshFlag 	    = True
        self.ApplyLineSearches      = True

        #definition of the convergence criteria
        #self.conv_criteria = DisplacementCriteria(0.000001,1e-9)

    #######################################################################
    def Initialize(self):

        self.time_scheme = ResidualBasedPredictorCorrectorBossakRotationScheme(self.damp_factor)

        #definition of the convergence criteria
        #self.conv_criteria = DisplacementCriteria(self.toll,self.absolute_tol)
        #builder_and_solver = ResidualBasedEliminationBuilderAndSolver(self.structure_linear_solver)

        #creating the solution strategy
	self.solver = ResidualBasedNewtonRaphsonLineSearchesStrategy(self.model_part,self.time_scheme,self.structure_linear_solver,self.conv_criteria,
	self.MaxNewtonRapshonIterations,self.MaxLineSearchIterations, self.tolls, self.amp, self.etmxa, self.etmna,
        self.CalculateReactionFlag,
        self.ReformDofSetAtEachStep,
        self.MoveMeshFlag,
        self.ApplyLineSearches)    
                 
    #######################################################################   
    def Solve(self):
        (self.solver).Solve()

    #######################################################################   
    def SetEchoLevel(self,level):
        (self.solver).SetEchoLevel(level)
