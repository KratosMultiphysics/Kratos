from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# importing the Kratos Library
from KratosMultiphysics import *
from KratosMultiphysics.StructuralApplication import *
from KratosMultiphysics.MixedElementApplication import *
from KratosMultiphysics.MeshingApplication import *
from KratosMultiphysics.ExternalSolversApplication import *


def AddVariables(model_part):
    model_part.AddNodalSolutionStepVariable(DISPLACEMENT)
    model_part.AddNodalSolutionStepVariable(ACCELERATION)
    model_part.AddNodalSolutionStepVariable(REACTION)
    model_part.AddNodalSolutionStepVariable(SX)
    model_part.AddNodalSolutionStepVariable(SY)
    model_part.AddNodalSolutionStepVariable(SZ)
    model_part.AddNodalSolutionStepVariable(SXY);
    model_part.AddNodalSolutionStepVariable(SXZ);
    model_part.AddNodalSolutionStepVariable(SYZ);
    model_part.AddNodalSolutionStepVariable(RADIATIVE_INTENSITY_1);  # radiative intensity is used as a commodity var!
    model_part.AddNodalSolutionStepVariable(RADIATIVE_INTENSITY_2);
    model_part.AddNodalSolutionStepVariable(RADIATIVE_INTENSITY_3);
    model_part.AddNodalSolutionStepVariable(RADIATIVE_INTENSITY_4);
    model_part.AddNodalSolutionStepVariable(RADIATIVE_INTENSITY_5);
    model_part.AddNodalSolutionStepVariable(RADIATIVE_INTENSITY_6);
    # model_part.AddNodalSolutionStepVariable(BODY_FORCE);
    model_part.AddNodalSolutionStepVariable(NEGATIVE_FACE_PRESSURE);
    model_part.AddNodalSolutionStepVariable(POSITIVE_FACE_PRESSURE);
    model_part.AddNodalSolutionStepVariable(DAMAGE);
    model_part.AddNodalSolutionStepVariable(NODAL_H);
    model_part.AddNodalSolutionStepVariable(FACE_LOAD);

    print("variables for the dynamic structural solution added correctly")


def AddDofs(model_part):
    for node in model_part.Nodes:
        # adding dofs
        node.AddDof(DISPLACEMENT_X, REACTION_X);
        node.AddDof(DISPLACEMENT_Y, REACTION_Y);
        node.AddDof(DISPLACEMENT_Z, REACTION_Z);
        node.AddDof(SX, RADIATIVE_INTENSITY_1)  # radiative intensity is used as a commodity var!
        node.AddDof(SY, RADIATIVE_INTENSITY_2)
        node.AddDof(SZ, RADIATIVE_INTENSITY_3)
        node.AddDof(SXY, RADIATIVE_INTENSITY_4)
        node.AddDof(SXZ, RADIATIVE_INTENSITY_5)
        node.AddDof(SYZ, RADIATIVE_INTENSITY_6)
    print("dofs for the dynamic structural solution added correctly")


class StaticStructuralSolver:
    #

    def __init__(self, model_part, domain_size):

        self.model_part = model_part
        self.time_scheme = ResidualBasedIncrementalUpdateStaticScheme()
        # if called here, Check may be called before the system is completely set up!!!
        self.time_scheme.Check(self.model_part)

        # definition of the solvers
        self.structure_linear_solver = SuperLUSolver()  # MKLPardisoSolver() #SkylineLUFactorizationSolver()

        # definition of the convergence criteria
#        self.conv_criteria = DisplacementCriteria(0.0001,1e-6)
        self.conv_criteria = MixedElementConvergenceCriteria(1e-2, 1e-20)
        self.conv_criteria.Check(self.model_part)
        self.MaxNewtonRapshonIterations = 100

        self.CalculateReactionFlag = True
        self.ReformDofSetAtEachStep = False
        self.MoveMeshFlag = False

        self.plot_util = PlotUtils()

        #(FindNodalHProcess(model_part)).Execute()

    #
    def Initialize(self):
        # creating the solution strategy

       # import strategy_python
       # self.solver = strategy_python.SolvingStrategyPython(self.model_part,self.time_scheme,self.structure_linear_solver,self.conv_criteria,self.CalculateReactionFlag,self.ReformDofSetAtEachStep,self.MoveMeshFlag)
        # (self.solver).SetEchoLevel(2)

        # creating the solution strategy
        self.solver = ResidualBasedNewtonRaphsonStrategy(self.model_part, self.time_scheme, self.structure_linear_solver, self.conv_criteria, self.MaxNewtonRapshonIterations, self.CalculateReactionFlag, self.ReformDofSetAtEachStep, self.MoveMeshFlag)
        self.solver.Check();

        #(self.solver).SetReformDofSetAtEachStepFlag(True)
        #(self.solver).SetMoveMeshFlag(True)

    #
    def Solve(self):
        (self.solver).Solve()

        self.plot_util.PlotVariable(DAMAGE, self.model_part)

        aux = Matrix(1, 6);
        aux2 = Matrix(1, 6);
        for node in self.model_part.Nodes:
            aux[0, 0] = node.GetSolutionStepValue(SX)
            aux[0, 1] = node.GetSolutionStepValue(SY)
            aux[0, 2] = node.GetSolutionStepValue(SZ)
            aux[0, 3] = 0.5 * node.GetSolutionStepValue(SXY)
            aux[0, 5] = 0.5 * node.GetSolutionStepValue(SXZ)  # note that the order of 4 and 5 is changed to match gid convenction
            aux[0, 4] = 0.5 * node.GetSolutionStepValue(SYZ)
            node.SetValue(GREEN_LAGRANGE_STRAIN_TENSOR, aux)

            d = node.GetSolutionStepValue(DAMAGE)

            for i in range(0, 6):
                aux2[0, i] = aux[0, i] * (1 - d);
            node.SetValue(CAUCHY_STRESS_TENSOR, aux2)

    #
    def SetEchoLevel(self, level):
        (self.solver).SetEchoLevel(level)
