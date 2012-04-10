#importing the Kratos Library
from KratosMultiphysics import *
from KratosMultiphysics.StructuralApplication import *
from KratosMultiphysics.TrilinosApplication import *
# Check that KratosMultiphysics was imported in the main script
CheckForPreviousImport()

def AddVariables(model_part):
    model_part.AddNodalSolutionStepVariable(DISPLACEMENT);
##    model_part.AddNodalSolutionStepVariable(DISPLACEMENT_OLD);
##    model_part.AddNodalSolutionStepVariable(DISPLACEMENT_NULL);
##    model_part.AddNodalSolutionStepVariable(DISPLACEMENT_EINS);
##    model_part.AddNodalSolutionStepVariable(DISPLACEMENT_DT);
##    model_part.AddNodalSolutionStepVariable(DISPLACEMENT_NULL_DT);
##    model_part.AddNodalSolutionStepVariable(DISPLACEMENT_EINS_DT);
##    model_part.AddNodalSolutionStepVariable(ACCELERATION);
##    model_part.AddNodalSolutionStepVariable(ACCELERATION_NULL);
##    model_part.AddNodalSolutionStepVariable(ACCELERATION_EINS);
##    model_part.AddNodalSolutionStepVariable(VELOCITY);
##    model_part.AddNodalSolutionStepVariable(ACCELERATION);
    model_part.AddNodalSolutionStepVariable(REACTION);
    model_part.AddNodalSolutionStepVariable(NEGATIVE_FACE_PRESSURE);
    model_part.AddNodalSolutionStepVariable(POSITIVE_FACE_PRESSURE);
    model_part.AddNodalSolutionStepVariable(INSITU_STRESS);
    model_part.AddNodalSolutionStepVariable(FACE_LOAD);
    model_part.AddNodalSolutionStepVariable(PARTITION_INDEX);

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
        self.time_scheme = TrilinosResidualBasedIncrementalUpdateStaticScheme()

        self.Comm = CreateCommunicator()

	self.buildertype="standard"

        #definition of the solvers
        self.structure_linear_solver =  TrilinosLinearSolver()
        
        #definition of the convergence criteria
        self.conv_criteria = TrilinosDisplacementCriteria(1e-6,1e-9,self.Comm)

        self.CalculateReactionFlag = False
        self.ReformDofSetAtEachStep = False
        self.MoveMeshFlag = True
        self.calculate_norm_dx_flag = False
        self.max_iterations = 10

        if(domain_size == 2):
            self.guess_row_size = 20
        else:
            self.guess_row_size = 45
            

        self.guess_row_size = 18
        
        
        
    #######################################################################
    def Initialize(self):
 
##        p_builder = TrilinosBuilderAndSolver(self.Comm,self.guess_row_size,self.structure_linear_solver)
## 
##        #creating the solution strategy
##        self.solver = ResidualBasedNewtonRaphsonStrategy(self.model_part, self.time_scheme,self.structure_linear_solver, self.conv_criteria,  p_builder,self.max_iterations, self.CalculateReactionFlag,self.ReformDofSetAtEachStep, self.MoveMeshFlag)
        import trilinos_strategy_python
        self.solver = trilinos_strategy_python.SolvingStrategyPython(self.buildertype,self.model_part,self.time_scheme,self.structure_linear_solver,self.conv_criteria,self.CalculateReactionFlag,self.ReformDofSetAtEachStep,self.MoveMeshFlag,self.Comm,self.guess_row_size)
        
 
                 
    #######################################################################   
    def Solve(self):
        (self.solver).Solve()

    #######################################################################   
    def SetEchoLevel(self,level):
        (self.solver).SetEchoLevel(level)

      

