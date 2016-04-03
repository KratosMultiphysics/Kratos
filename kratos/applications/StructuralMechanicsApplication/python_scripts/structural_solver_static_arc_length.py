from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# importing the Kratos Library
from KratosMultiphysics import *
from KratosMultiphysics.SolidMechanicsApplication import *
from KratosMultiphysics.StructuralMechanicsApplication import *
CheckForPreviousImport()

def AddVariables(model_part, config=None):
    model_part.AddNodalSolutionStepVariable(NODAL_MASS) #MSI, i included the variable becouse i calculate energy
    # add displacements
    model_part.AddNodalSolutionStepVariable(DISPLACEMENT)
    # add dynamic variables
    model_part.AddNodalSolutionStepVariable(VELOCITY)
    model_part.AddNodalSolutionStepVariable(ACCELERATION)
    # add reactions for the displacements
    model_part.AddNodalSolutionStepVariable(REACTION)
    # add nodal force variables
    model_part.AddNodalSolutionStepVariable(INTERNAL_FORCE)
    model_part.AddNodalSolutionStepVariable(EXTERNAL_FORCE)
    model_part.AddNodalSolutionStepVariable(CONTACT_FORCE)
    # add specific variables for the problem conditions
    model_part.AddNodalSolutionStepVariable(POSITIVE_FACE_PRESSURE)
    model_part.AddNodalSolutionStepVariable(NEGATIVE_FACE_PRESSURE)
    model_part.AddNodalSolutionStepVariable(POINT_LOAD)
    model_part.AddNodalSolutionStepVariable(LINE_LOAD)
    model_part.AddNodalSolutionStepVariable(SURFACE_LOAD)
    model_part.AddNodalSolutionStepVariable(VOLUME_ACCELERATION)

    print("*********************************************************************** ")
    print("Variables for the Static Structural Arc Length Solution added correctly")
    print("*********************************************************************** ")


def ChangeCondition(model_part, lamda):
    for node in model_part.Nodes:
        new_load = node.GetSolutionStepValue(POINT_LOAD) * lamda;
        node.SetSolutionStepValue(POINT_LOAD, 0, new_load)

def AddDofs(model_part, config=None):
    for node in model_part.Nodes:
    # adding dofs
        node.AddDof(DISPLACEMENT_X, REACTION_X);
        node.AddDof(DISPLACEMENT_Y, REACTION_Y);
        node.AddDof(DISPLACEMENT_Z, REACTION_Z);

    print("*********************************************************************** ")
    print("Dofs for the Static Structural Arc Length Solution added correctly")
    print("*********************************************************************** ")


class StaticArcLengthStructuralSolver:
    #
    def __init__(self, model_part, domain_size):

        self.model_part = model_part
        self.time_scheme = ResidualBasedIncrementalUpdateStaticScheme()

        # self.time_scheme   = ParallelResidualBasedIncrementalUpdateStaticScheme()
        # Varibles de Control de Arc Lenght Method
        self.Ide = 5
        self.factor_delta_lmax = 1.00
        self.toler = 1.0E-9
        self.norm = 1.0E-5
        self.max_iter = 20

         # self.time_scheme = ResidualBasedPredictorCorrectorBossakScheme(self.damp_factor)
        # definition of the solvers. Super_Lu Default
        self.structure_linear_solver = SkylineLUFactorizationSolver()
        # self.structure_linear_solver      =   SuperLUSolver()
        # self.structure_linear_solver =  MKLPardisoSolver()
        # pDiagPrecond = ParallelDiagonalPreconditioner()
        # self.structure_linear_solver     =  ParallelCGSolver(1e-8, 5000,pDiagPrecond)
        # self.structure_linear_solver    =   Preconditioner()
        # self.structure_linear_solver    =   IterativeSolver()

        # pDiagPrecond = DiagonalPreconditioner()
        # LST  = 1E-9
        # LSMI = 5000
        # self.structure_linear_solver  =  BICGSTABSolver(LST,LSMI,pDiagPrecond)

        # definition of the convergence criteria
        Displacement = DisplacementCriteria(1E-6, 1E-9)
        Residual = ResidualCriteria(1E-3, 1E-6)

        self.conv_criteria = AndCriteria(Residual, Displacement)
        # self.conv_Residual     = ResidualCriteria(0.000001,1E-9)
        # self.conv_Displacement = DisplacementCriteria(self.norm,self.toler)
        # self.conv_criteria     = ResDisCriteria(self.conv_Residual, self.conv_Displacement)
        # self.conv_criteria     = ResidualDisplacementCriteria(self.norm,self.toler)
        # self.conv_criteria     = DisplacementCriteria(self.norm,self.toler)
        # self.conv_criteria    = ParallelDisplacementCriteria(0.000001,1e-9)
        # self.conv_criteria     = DisplacementCriteria(self.norm,self.toler)

        # definition of the convergence criteria

        self.CalculateReactionFlag = True
        self.ReformDofSetAtEachStep = True
        self.MoveMeshFlag = True

    #
    def Initialize(self):

        self.solver = ResidualBasedArcLengthStrategy(
            self.model_part, 
            self.time_scheme,
            self.structure_linear_solver, 
            self.conv_criteria, 
            self.Ide, 
            self.max_iter, 
            self.factor_delta_lmax, 
            self.CalculateReactionFlag, 
            self.ReformDofSetAtEachStep, 
            self.MoveMeshFlag
            )
    #
    def Solve(self):
        (self.solver).Solve()

    #
    def SetEchoLevel(self, level):
        (self.solver).SetEchoLevel(level)

    def ChangeCondition(self, model_part, lamda):
        for cond in model_part.Conditions:
            print(cond)
            
def CreateSolver(model_part, config):
    structural_solver = StaticArcLengthStructuralSolver(model_part, config.domain_size)
    model_part.ProcessInfo[LAMBDA] = 0.00;
    return structural_solver
