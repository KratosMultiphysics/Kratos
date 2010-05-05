#importing the Kratos Library
from Kratos import *
from KratosIncompressibleFluidApplication import *
from KratosTrilinosApplication import *


def AddVariables(model_part):
    model_part.AddNodalSolutionStepVariable(VELOCITY);
    model_part.AddNodalSolutionStepVariable(MESH_VELOCITY);
    model_part.AddNodalSolutionStepVariable(PRESSURE);
    model_part.AddNodalSolutionStepVariable(IS_FLUID);
    model_part.AddNodalSolutionStepVariable(IS_STRUCTURE);
    model_part.AddNodalSolutionStepVariable(IS_FREE_SURFACE);
    model_part.AddNodalSolutionStepVariable(IS_BOUNDARY);
    model_part.AddNodalSolutionStepVariable(DISPLACEMENT);
    model_part.AddNodalSolutionStepVariable(VISCOSITY);
    model_part.AddNodalSolutionStepVariable(DENSITY);
    model_part.AddNodalSolutionStepVariable(BODY_FORCE);
    model_part.AddNodalSolutionStepVariable(NODAL_AREA);
    model_part.AddNodalSolutionStepVariable(ADVPROJ)
    model_part.AddNodalSolutionStepVariable(DIVPROJ)
    model_part.AddNodalSolutionStepVariable(THAWONE)
    model_part.AddNodalSolutionStepVariable(THAWTWO)
    model_part.AddNodalSolutionStepVariable(REACTION); 
    model_part.AddNodalSolutionStepVariable(REACTION_WATER_PRESSURE); 
    print "variables for the dynamic structural solution added correctly"
        
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
##        self.time_scheme = ResidualBasedIncrementalUpdateStaticScheme()

        self.Comm = CreateCommunicator()

        self.move_mesh_strategy = 0
        self.time_scheme = TrilinosResidualBasedLagrangianMonolithicScheme( self.move_mesh_strategy)

        #definition of the solvers
        self.linear_solver =  TrilinosLinearSolver()
        
        #definition of the convergence criteria
##        self.conv_criteria = UPCriteria(1e-7,1e-9,1e-7,1e-9)
        self.conv_criteria = TrilinosUPCriteria(1e-3,1e-9,1e-3,1e-6,self.Comm)

        self.model_part.ProcessInfo.SetValue(DYNAMIC_TAU, 0.001);
        self.max_iter = 10
                            
        #default settings
        self.echo_level = 1
        self.CalculateReactionFlag = False
        self.ReformDofSetAtEachStep = False
        self.CalculateNormDxFlag = True
        self.MoveMeshFlag = False

        if(domain_size == 2):
            self.guess_row_size = 15
        else:
            self.guess_row_size = 45
            

        self.guess_row_size = 18

        self.buildertype="standard"
        
    #######################################################################
    def Initialize(self):
        #creating the solution strategy
        
        import trilinos_strategy_python
        self.solver = trilinos_strategy_python.SolvingStrategyPython(self.buildertype,self.model_part,self.time_scheme,self.linear_solver,self.conv_criteria,self.CalculateReactionFlag,self.ReformDofSetAtEachStep,self.MoveMeshFlag,self.Comm,self.guess_row_size)


##        (self.neigh_finder).Execute();
##        self.Remesh()        
                 
    #######################################################################   
    def Solve(self):
        (self.solver).Solve()


        

    #######################################################################   
    def SetEchoLevel(self,level):
        (self.solver).SetEchoLevel(level)
    
    ########################################################################





