from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# importing the Kratos Library
from KratosMultiphysics import *
from KratosMultiphysics.StructuralApplication import *
CheckForPreviousImport()


def AddVariables(model_part):
    model_part.AddNodalSolutionStepVariable(FORCE)
    model_part.AddNodalSolutionStepVariable(DISPLACEMENT)
    model_part.AddNodalSolutionStepVariable(NODAL_MASS)
    model_part.AddNodalSolutionStepVariable(RHS)
    model_part.AddNodalSolutionStepVariable(VELOCITY)
    model_part.AddNodalSolutionStepVariable(ACCELERATION)
    model_part.AddNodalSolutionStepVariable(REACTION)
    model_part.AddNodalSolutionStepVariable(NEGATIVE_FACE_PRESSURE)
    model_part.AddNodalSolutionStepVariable(POSITIVE_FACE_PRESSURE)
    # model_part.AddNodalSolutionStepVariable(BODY_FORCE);
    model_part.AddNodalSolutionStepVariable(DAMAGE);
    # model_part.AddNodalSolutionStepVariable(DENSITY);
    # model_part.AddNodalSolutionStepVariable(PRESSURE);
    model_part.AddNodalSolutionStepVariable(DISPLACEMENT_OLD);
    model_part.AddNodalSolutionStepVariable(NODAL_VALUES);
    model_part.AddNodalSolutionStepVariable(IS_BOUNDARY);
    model_part.AddNodalSolutionStepVariable(IS_CONTACT_SLAVE);
    model_part.AddNodalSolutionStepVariable(IS_CONTACT_MASTER);
    model_part.AddNodalSolutionStepVariable(NODAL_AREA);
    model_part.AddNodalSolutionStepVariable(NODAL_VOLUME);
    model_part.AddNodalSolutionStepVariable(NODAL_DAMAGE);
    model_part.AddNodalSolutionStepVariable(IS_DUPLICATED);
    model_part.AddNodalSolutionStepVariable(NODAL_STRAIN);
    model_part.AddNodalSolutionStepVariable(REFINEMENT_LEVEL);
    model_part.AddNodalSolutionStepVariable(SPLIT_NODAL);
    model_part.AddNodalSolutionStepVariable(NORMAL);
    model_part.AddNodalSolutionStepVariable(GRAVITY);
    model_part.AddNodalSolutionStepVariable(POTENCIAL_ENERGY);
    model_part.AddNodalSolutionStepVariable(KINETIC_ENERGY);
    model_part.AddNodalSolutionStepVariable(INTERFACE_FORCES);
    model_part.AddNodalSolutionStepVariable(NORMAL_CONTACT_STRESS);
    model_part.AddNodalSolutionStepVariable(GRAVITY);
    model_part.AddNodalSolutionStepVariable(POTENCIAL_ENERGY);
    model_part.AddNodalSolutionStepVariable(KINETIC_ENERGY);
    model_part.AddNodalSolutionStepVariable(NODAL_STRESS);
    print("Variables for the dynamic structural solution added correctly")


def AddDofs(model_part):
    for node in model_part.Nodes:
        node.AddDof(DISPLACEMENT_X, REACTION_X);
        node.AddDof(DISPLACEMENT_Y, REACTION_Y);
        node.AddDof(DISPLACEMENT_Z, REACTION_Z);
    print("Dofs for the dynamic structural solution added correctly")


class DynamicStructuralSolver:
    #

    def __init__(self, model_part, domain_size):

        self.model_part = model_part
        self.domain_size = domain_size
        self.damping_ratio = 0.00;
        self.penalty_factor = 10.00
        self.max_delta_time = 0.05;
        self.fraction_delta_time = 0.90;
        self.CalculateReactionFlag = True;
        self.MoveMeshFlag = True;
        self.ComputeContactConditions = False;
        self.CE = Constraint_Enforcement.Penalty_Methods;  # Constraint_Enforcement.Lagrange_Multiplie_Methods;
        self.structure_linear_solver = SkylineLUFactorizationSolver()
        self.time_scheme = ResidualBasedIncrementalUpdateStaticScheme()
        self.builder = ResidualBasedEliminationBuilderAndSolver(self.structure_linear_solver);
        # self.time_scheme.Check(self.model_part)
        # self.structure_linear_solver =   SuperLUSolver()
        # self.structure_linear_solver  =  MKLPardisoSolver()
        # pDiagPrecond = ParallelDiagonalPreconditioner()
        # self.structure_linear_solver    =  ParallelCGSolver(1e-8, 5000,pDiagPrecond)
        # self.structure_linear_solver    =   Preconditioner()
        # self.structure_linear_solver    =   IterativeSolver()

    #

    def CriticalTime(self):
        (self.solver).Initialize();
        print("Calculating Time Step ")
        (self.solver).ComputeCriticalTime()

    def Initialize(self):

        # creating the solution strategy
        self.solver = ResidualBasedCentralDiferencesStrategy(self.model_part, self.CE, self.domain_size, self.damping_ratio, self.fraction_delta_time, self.max_delta_time, self.penalty_factor, self.CalculateReactionFlag, self.ComputeContactConditions, self.MoveMeshFlag, self.structure_linear_solver, self.time_scheme, self.builder)

    #
    def Solve(self):
        (self.solver).Solve()

    #
    def SetEchoLevel(self, level):
        (self.solver).SetEchoLevel(level)

    def SetFractionDeltaTime(self, fraction):
        (self.solver).SetFractionDeltaTime(fraction)

    def SetConditionsFlag(self, ComputeContactConditions):
        (self.solver).SetConditionsFlag(ComputeContactConditions)

    def CalculateBoundaryContours(self, ComputeBoundary):
        (self.solver).CalculateBoundaryContours(ComputeBoundary)
