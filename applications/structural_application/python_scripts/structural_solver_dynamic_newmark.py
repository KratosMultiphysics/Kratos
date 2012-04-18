#importing the Kratos Library
from KratosMultiphysics import *
from KratosMultiphysics.StructuralApplication import *
CheckForPreviousImport()



def AddVariables(model_part):
    model_part.AddNodalSolutionStepVariable(DISPLACEMENT);
    model_part.AddNodalSolutionStepVariable(VELOCITY);
    model_part.AddNodalSolutionStepVariable(ACCELERATION);
    model_part.AddNodalSolutionStepVariable(REACTION);
    model_part.AddNodalSolutionStepVariable(NEGATIVE_FACE_PRESSURE);
    model_part.AddNodalSolutionStepVariable(POSITIVE_FACE_PRESSURE);
    model_part.AddNodalSolutionStepVariable(DISPLACEMENT_OLD);
    model_part.AddNodalSolutionStepVariable(DISPLACEMENT_NULL);
    model_part.AddNodalSolutionStepVariable(DISPLACEMENT_EINS);
    model_part.AddNodalSolutionStepVariable(DISPLACEMENT_DT);
    model_part.AddNodalSolutionStepVariable(DISPLACEMENT_NULL_DT);
    model_part.AddNodalSolutionStepVariable(DISPLACEMENT_EINS_DT);
    model_part.AddNodalSolutionStepVariable(ACCELERATION_NULL);
    model_part.AddNodalSolutionStepVariable(ACCELERATION_EINS);
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
        
        self.damp_factor = 1.0
        self.toll = 1e-6
        self.absolute_tol = 1e-9

        #definition of the solvers
        self.structure_linear_solver =  SkylineLUFactorizationSolver()
        
        #definition of the convergence criteria
        #self.conv_criteria = DisplacementCriteria(0.000001,1e-9)

    #######################################################################
    def Initialize(self):

        self.time_scheme = ResidualBasedNewmarkScheme(self.damp_factor)

        #definition of the convergence criteria
        self.conv_criteria = DisplacementCriteria(self.toll,self.absolute_tol)
    
        #creating the solution strategy
        self.solver = ResidualBasedNewtonRaphsonStrategy(self.model_part,self.time_scheme,self.structure_linear_solver,self.conv_criteria,30,True,True,True)
        print "self.echo_level = " , self.echo_level
        (self.solver).SetEchoLevel(self.echo_level)
        #(self.solver).SetMoveMeshFlag(True)
        #(self.solver).SetReformDofSetAtEachStepFlag(True)
        #(self.solver).SetBuilderAndSolver(builder_and_solver)
        print "finished initialization of the fluid strategy"

                 
    #######################################################################   
    def Solve(self):
        (self.solver).Solve()

