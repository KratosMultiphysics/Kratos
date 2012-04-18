#importing the Kratos Library
from KratosMultiphysics import *
from KratosMultiphysics.StructuralApplication import *
from KratosMultiphysics.ExternalSolversApplication import *
CheckForPreviousImport()


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

class StaticStructuralSolver:
    #######################################################################
    def __init__(self,model_part,domain_size):

        self.model_part = model_part
        self.time_scheme = ResidualBasedIncrementalUpdateStaticScheme()

        #definition of the solvers
        #self.structure_linear_solver =  SkylineLUFactorizationSolver()
        #preconditioner = ILU0Preconditioner()
        preconditioner = DiagonalPreconditioner()
        #self.structure_linear_solver = BICGSTABSolver(1e-12, 15000, preconditioner)
        #self.structure_linear_solver = MGMRESSolver(1e-12, 15000, preconditioner)
        #self.structure_linear_solver = TFQMRSolver(1e-8, 15000, preconditioner)
        self.structure_linear_solver = MGMRESSolver(1e-6, 15000 )
        #definition of the convergence criteria
        self.conv_criteria = DisplacementCriteria(0.000001,1e-9)

    #######################################################################
    def Initialize(self):
        #creating the solution strategy
        
        #builder_and_solver = ResidualBasedEliminationBuilderAndSolver(self.structure_linear_solver)
        #builder_and_solver = DeactivatingBuilderAndSolver(self.structure_linear_solver)
        
        self.solver = ResidualBasedUzawaNewtonRaphsonStrategy(self.model_part,self.time_scheme,self.structure_linear_solver,self.conv_criteria,30,True,True,True)
        #(self.solver).SetMoveMeshFlag(True)
        #(self.solver).SetReformDofSetAtEachStepFlag(True)
        #(self.solver).SetBuilderAndSolver(builder_and_solver)
        
 
                 
    #######################################################################   
    def Solve(self):
        (self.solver).Solve()

