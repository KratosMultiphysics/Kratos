from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# importing the Kratos Library
from KratosMultiphysics import *
from KratosMultiphysics.ShallowWaterApplication import *
CheckForPreviousImport()


def AddVariables(model_part):  #this way er only need one command to add all the variables to our problem 
    model_part.AddNodalSolutionStepVariable(VELOCITY);
    model_part.AddNodalSolutionStepVariable(HEIGHT);
    # Physic problem parameters
    model_part.AddNodalSolutionStepVariable(BATHYMETRY);
    model_part.AddNodalSolutionStepVariable(GRAVITY);
    model_part.AddNodalSolutionStepVariable(RAIN);

def AddDofs(model_part):
    for node in model_part.Nodes:
        node.AddDof(VELOCITY_X);
        node.AddDof(VELOCITY_Y);
        node.AddDof(HEIGHT);
    print ("variables for the SWE solver added correctly")

class ShallowWaterSolver:
    #######################################################################
    def __init__(self,model_part,domain_size):  # Constructor of the class
        self.model_part = model_part
        self.domain_size = domain_size
        self.time_scheme = ResidualBasedIncrementalUpdateStaticScheme()

        # Definition of the mesh stage solver
        self.swe_linear_solver =  SkylineLUFactorizationSolver()  # We set the type of solver that we want
        
        # definition of the convergence criteria
        self.conv_criteria = DisplacementCriteria(1e-6,1e-9)  # Tolerance for the solver
    
    #######################################################################
    def Initialize(self):
        # Creating the solution strategy for the mesh stage
        CalculateReactionFlag = False
        ReformDofSetAtEachStep = False
        MoveMeshFlag = False
        import strategy_python
        self.mesh_solver = strategy_python.SolvingStrategyPython(self.model_part,
                                            self.time_scheme, self.swe_linear_solver,
                                            self.conv_criteria, CalculateReactionFlag,
                                            ReformDofSetAtEachStep, MoveMeshFlag)
        # Initialize dry/wet state utility
        self.drybedutility = DryBedUtility(self.model_part)

    #######################################################################
    def Solve(self):        
        # Check dry bed
        (self.drybedutility).CheckPrimitiveVariables(self.model_part.Nodes)
        # Solve equations on mesh
        (self.mesh_solver).Solve()
        
    #######################################################################
    def SetEchoLevel(self,level):
        (self.mesh_solver).SetEchoLevel(level)
