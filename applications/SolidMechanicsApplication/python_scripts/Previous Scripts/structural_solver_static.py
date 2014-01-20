from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# importing the Kratos Library
from KratosMultiphysics import *
from KratosMultiphysics.SolidMechanicsApplication import *
CheckForPreviousImport()


def AddVariables(model_part):
    model_part.AddNodalSolutionStepVariable(DISPLACEMENT)
    model_part.AddNodalSolutionStepVariable(PRESSURE)
    model_part.AddNodalSolutionStepVariable(ACCELERATION)
    model_part.AddNodalSolutionStepVariable(VELOCITY)
    model_part.AddNodalSolutionStepVariable(ACCELERATION)
    model_part.AddNodalSolutionStepVariable(REACTION)
    model_part.AddNodalSolutionStepVariable(FORCE_INTERNAL)
    model_part.AddNodalSolutionStepVariable(FORCE_EXTERNAL)
    model_part.AddNodalSolutionStepVariable(FORCE_DYNAMIC)
    model_part.AddNodalSolutionStepVariable(NEGATIVE_FACE_PRESSURE);
    model_part.AddNodalSolutionStepVariable(POSITIVE_FACE_PRESSURE);
    model_part.AddNodalSolutionStepVariable(FACE_LOAD);

    print("variables for the dynamic structural solution added correctly")


def AddDofs(model_part):
    for node in model_part.Nodes:
        # adding dofs
        node.AddDof(DISPLACEMENT_X, REACTION_X);
        node.AddDof(DISPLACEMENT_Y, REACTION_Y);
        node.AddDof(DISPLACEMENT_Z, REACTION_Z);
        node.AddDof(PRESSURE, REACTION_PRESSURE);
    print("dofs for the dynamic structural solution added correctly")


class StaticStructuralSolver:
    #

    def __init__(self, model_part, domain_size):

        self.model_part = model_part
        self.time_scheme = ResidualBasedStaticScheme()
        # if called here, Check may be called before the system is completely set up!!!
        self.time_scheme.Check(self.model_part)

        # definition of the solvers
        self.structure_linear_solver = SkylineLUFactorizationSolver()

        # definition of the convergence criteria
#        self.conv_criteria = DisplacementCriteria(0.0001,1e-6)
        self.conv_criteria = DisplacementCriteria(1e-6, 1e-9)
        self.conv_criteria.Check(self.model_part)
        self.MaxNewtonRapshonIterations = 100

        self.CalculateReactionFlag = True
        self.ReformDofSetAtEachStep = False
        self.MoveMeshFlag = True

    #
    def Initialize(self):
        # creating the solution strategy

       # import strategy_python
       # self.solver = strategy_python.SolvingStrategyPython(self.model_part,self.time_scheme,self.structure_linear_solver,self.conv_criteria,self.CalculateReactionFlag,self.ReformDofSetAtEachStep,self.MoveMeshFlag)
        # (self.solver).SetEchoLevel(2)

        # creating the solution strategy
        self.solver = ResidualBasedNewtonRaphsonStrategy(self.model_part, self.time_scheme, self.structure_linear_solver, self.conv_criteria, self.MaxNewtonRapshonIterations, self.CalculateReactionFlag, self.ReformDofSetAtEachStep, self.MoveMeshFlag)

        # solving with AMG
        # builder_and_solver = BlockResidualBasedBuilderAndSolver(self.structure_linear_solver)
        # self.solver = ResidualBasedNewtonRaphsonStrategy(self.model_part,self.time_scheme,self.structure_linear_solver,self.conv_criteria,builder_and_solver,self.MaxNewtonRapshonIterations,self.CalculateReactionFlag,self.ReformDofSetAtEachStep,self.MoveMeshFlag)

        self.solver.Check();
        #(self.solver).SetReformDofSetAtEachStepFlag(True)
        #(self.solver).SetMoveMeshFlag(True)

    #
    def Solve(self):
        (self.solver).Solve()

    #
    def SetEchoLevel(self, level):
        (self.solver).SetEchoLevel(level)
