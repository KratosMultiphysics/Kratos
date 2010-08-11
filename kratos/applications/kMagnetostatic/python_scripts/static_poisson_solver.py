#importing the Kratos Library
from Kratos import *
from KratosR1MagnetostaticApplication import *

def AddVariables(model_part):
    model_part.AddNodalSolutionStepVariable(DISPLACEMENT);
    model_part.AddNodalSolutionStepVariable(MAGNETIC_PERMEABILITY);
    model_part.AddNodalSolutionStepVariable(COERCIVITY);
    model_part.AddNodalSolutionStepVariable(MAGNETIC_FIELD_INTENSITY);
    model_part.AddNodalSolutionStepVariable(MAGNETIC_FLUX_DENSITY);

    model_part.AddNodalSolutionStepVariable(MAGNETOSTATIC_POTENTIAL);
    model_part.AddNodalSolutionStepVariable(MAGNETOSTATIC_VECTOR_POTENTIAL);
    model_part.AddNodalSolutionStepVariable(MAGNETOSTATIC_POINT_CURRENT);
    model_part.AddNodalSolutionStepVariable(INFINIT_COEFFICIENT);

def AddDofs(model_part):
    for node in model_part.Nodes:

        #adding dofs
        node.AddDof(MAGNETOSTATIC_POTENTIAL);
        node.AddDof(MAGNETOSTATIC_VECTOR_POTENTIAL_X);
        node.AddDof(MAGNETOSTATIC_VECTOR_POTENTIAL_Y);
        node.AddDof(MAGNETOSTATIC_VECTOR_POTENTIAL_Z);

    print "variables for the Magnetostatic solver added correctly"

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
        MoveMeshFlag = False
        import strategy_python
        self.solver = strategy_python.SolvingStrategyPython(self.model_part,self.time_scheme,self.poisson_linear_solver,self.conv_criteria,CalculateReactionFlag,ReformDofSetAtEachStep,MoveMeshFlag)
      
                 
    #######################################################################   
    def Solve(self):
        (self.solver).Solve()

    #######################################################################   
    def SetEchoLevel(self,level):
        (self.solver).SetEchoLevel(level)
