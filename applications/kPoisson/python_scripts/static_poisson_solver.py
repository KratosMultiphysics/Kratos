#importing the Kratos Library
from Kratos import *
from KratosR1PoissonApplication import *

def AddVariables(model_part):
    model_part.AddNodalSolutionStepVariable(DUMMY_UNKNOWN);
    model_part.AddNodalSolutionStepVariable(DUMMY_POINT_SOURCE);

def AddDofs(model_part):
    for node in model_part.Nodes:

        #adding dofs
        node.AddDof(DUMMY_UNKNOWN);

    print "variables for the Poisson solver added correctly"

class StaticPoissonSolver:
    #######################################################################
    def __init__(self,model_part,domain_size):

        self.model_part = model_part
        self.time_scheme = ResidualBasedIncrementalUpdateStaticScheme()

        #definition of the solvers
        self.poisson_linear_solver =  SkylineLUFactorizationSolver()
        
        #definition of the convergence criteria
        self.conv_criteria = DisplacementCriteria(1e-6,1e-9)
        
    #######################################################################
    def Initialize(self):
        #creating the solution strategy
        CalculateReactionFlag = False
        ReformDofSetAtEachStep = False
        MoveMeshFlag = True
        import strategy_python
        self.solver = strategy_python.SolvingStrategyPython(self.model_part,self.time_scheme,self.poisson_linear_solver,self.conv_criteria,CalculateReactionFlag,ReformDofSetAtEachStep,MoveMeshFlag)
      
                 
    #######################################################################   
    def Solve(self):
        (self.solver).Solve()

    #######################################################################   
    def SetEchoLevel(self,level):
        (self.solver).SetEchoLevel(level)
