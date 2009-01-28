#importing the Kratos Library
from Kratos import *
from KratosStructuralApplication import *
from KratosExternalSolversApplication import *

def AddVariables(model_part):
    model_part.AddNodalSolutionStepVariable(DISPLACEMENT);
    model_part.AddNodalSolutionStepVariable(VELOCITY);
    model_part.AddNodalSolutionStepVariable(ACCELERATION);
    model_part.AddNodalSolutionStepVariable(REACTION);
    model_part.AddNodalSolutionStepVariable(NEGATIVE_FACE_PRESSURE);
    model_part.AddNodalSolutionStepVariable(POSITIVE_FACE_PRESSURE);
    model_part.AddNodalSolutionStepVariable(FACE_LOAD);

    print "variables for the dynamic structural solution added correctly"
        
def AddDofs(model_part):
    for node in model_part.Nodes:
        #adding dofs
        node.AddDof(DISPLACEMENT_X);
        node.AddDof(DISPLACEMENT_Y);
        node.AddDof(DISPLACEMENT_Z);
    print "dofs for the dynamic structural solution added correctly"


class DynamicStructuralSolver:
    #######################################################################
    def __init__(self,model_part,domain_size):

        self.model_part = model_part
        self.echo_level = 0
        
        self.damp_factor = -0.1
        self.toll = 1e-6
        self.absolute_tol = 1e-9

        #definition of the solvers
        self.structure_linear_solver =  SuperLUSolver()

        #definition of the convergence criteria
        self.conv_criteria = DisplacementCriteria(0.0001,1e-6)

    #######################################################################
    def Initialize(self):

        self.time_scheme = ResidualBasedPredictorCorrectorBossakScheme(self.damp_factor)

        #definition of the convergence criteria
        self.conv_criteria = DisplacementCriteria(self.toll,self.absolute_tol)
    
        #creating the solution strategy
        self.solver = ResidualBasedNewtonRaphsonStrategy(self.model_part,self.time_scheme,self.structure_linear_solver,self.conv_criteria)
#        print "self.echo_level = " , self.echo_level
#        (self.solver).SetEchoLevel(self.echo_level)
        (self.solver).SetMoveMeshFlag(True)
        print "finished initialization of the structural strategy"

                 
    #######################################################################   
    def Solve(self):
        (self.solver).Solve()

