#importing the Kratos Library
from KratosMultiphysics import *
from KratosMultiphysics.StructuralApplication import *
from KratosMultiphysics.MKLSolversApplication import *

def AddVariables(model_part):
    model_part.AddNodalSolutionStepVariable(DISPLACEMENT);
    model_part.AddNodalSolutionStepVariable(DISPLACEMENT_OLD);
    model_part.AddNodalSolutionStepVariable(DISPLACEMENT_NULL);
    model_part.AddNodalSolutionStepVariable(DISPLACEMENT_EINS);
    model_part.AddNodalSolutionStepVariable(DISPLACEMENT_DT);
    model_part.AddNodalSolutionStepVariable(DISPLACEMENT_NULL_DT);
    model_part.AddNodalSolutionStepVariable(DISPLACEMENT_EINS_DT);
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
        node.AddDof(DISPLACEMENT_X,REACTION_X);
        node.AddDof(DISPLACEMENT_Y,REACTION_Y);
        node.AddDof(DISPLACEMENT_Z,REACTION_Z);
    print "dofs for the dynamic structural solution added correctly"


class StaticStructuralSolver:
    #######################################################################
    def __init__(self,model_part,domain_size):

        self.model_part = model_part
        self.echo_level = 0
        
        self.damp_factor = -0.1
        self.toll = 1e-6
        self.absolute_tol = 1e-9

        self.CalculateReactionFlag = True
        self.ReformDofSetAtEachStep = False
        self.MoveMeshFlag = True

        #definition of the solvers
##        self.structure_linear_solver =  ParallelSkylineLUFactorizationSolver()
        pDiagPrecond = ParallelDiagonalPreconditioner()
        self.structure_linear_solver =  ParallelCGSolver(1e-8, 5000,pDiagPrecond)

        #definition of the convergence criteria
        self.conv_criteria = ParallelDisplacementCriteria(0.000001,1e-9)

    #######################################################################
    def Initialize(self):

        self.time_scheme = ParallelResidualBasedIncrementalUpdateStaticScheme()

        #definition of the convergence criteria
        self.conv_criteria = ParallelDisplacementCriteria(self.toll,self.absolute_tol)
        #builder_and_solver = ResidualBasedEliminationBuilderAndSolver(self.structure_linear_solver)
    
        #creating the solution strategy

        import parallel_strategy_python
        self.solver = parallel_strategy_python.SolvingStrategyPython(self.model_part,self.time_scheme,self.structure_linear_solver,self.conv_criteria,self.CalculateReactionFlag,self.ReformDofSetAtEachStep,self.MoveMeshFlag)
        
        
##        self.solver = ResidualBasedNewtonRaphsonStrategy(self.model_part,self.time_scheme,self.structure_linear_solver,self.conv_criteria,30,True,False,True)
        print "self.echo_level = " , self.echo_level
        #(self.solver).SetBuilderAndSolver(builder_and_solver)
##        (self.solver).SetEchoLevel(self.echo_level)
        #(self.solver).SetReformDofSetAtEachStepFlag(True)
        #(self.solver).SetMoveMeshFlag(True)
        print "finished initialization of the fluid strategy"

                 
    #######################################################################   
    def Solve(self):
        (self.solver).Solve()

    #######################################################################   
    def SetEchoLevel(self,level):
        (self.solver).SetEchoLevel(level)
