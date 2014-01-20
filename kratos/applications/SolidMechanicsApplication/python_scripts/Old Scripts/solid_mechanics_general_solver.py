from __future__ import unicode_literals, print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# importing the Kratos Library
from KratosMultiphysics import *
from KratosMultiphysics.SolidMechanicsApplication import *
CheckForPreviousImport()


def AddVariables(model_part):
    model_part.AddNodalSolutionStepVariable(DISPLACEMENT)
    model_part.AddNodalSolutionStepVariable(PRESSURE)
    model_part.AddNodalSolutionStepVariable(VELOCITY)
    model_part.AddNodalSolutionStepVariable(ACCELERATION)
    model_part.AddNodalSolutionStepVariable(REACTION)
    model_part.AddNodalSolutionStepVariable(FORCE_INTERNAL)
    model_part.AddNodalSolutionStepVariable(FORCE_EXTERNAL)
    model_part.AddNodalSolutionStepVariable(FORCE_DYNAMIC)
    model_part.AddNodalSolutionStepVariable(REACTION_PRESSURE)

    print("VARIABLES ADDED CORRECTLY")


def AddDofs(model_part, problemtype):
    if(problemtype == "Mechanical"):
        for node in model_part.Nodes:
            node.AddDof(DISPLACEMENT_X, REACTION_X)
            node.AddDof(DISPLACEMENT_Y, REACTION_Y);
            node.AddDof(DISPLACEMENT_Z, REACTION_Z);
            node.AddDof(PRESSURE, REACTION_PRESSURE);

    print("DOF'S ADDED CORRECTLY")


class SolidSolution:
    #

    def __init__(self, model_part, domain_size, problemtype, dynamicstype):

        self.model_part = model_part
        self.echo_level = 1

        # type of problem (mechanical, thermal, thermo-mechanical)
        self.problemtype = problemtype

        # type of time evolution (dynamic=1, static=0)
        self.dynamicstype = dynamicstype

        # definition of the solvers
        self.solid_linear_solver = SkylineLUFactorizationSolver()
        # self.solid_linear_solver      =   SuperLUSolver()

        # definition of the convergence criteria
        self.rel_tol = 1e-6
        self.abs_tol = 1e-9
        self.max_iter = 30;

        # mechanics computation
        self.ComputeMechanicsFlag = False;

        # self.mechanical_convergence_criteria = DisplacementConvergenceCriteria(self.abs_tol,self.abs_tol)
        self.mechanical_convergence_criteria = DisplacementCriteria(self.abs_tol, self.abs_tol)

        if(self.problemtype == "Mechanical" or self.problemtype == "ThermoMechanical"):
            self.mechanical_convergence_criteria.Check(self.model_part)

        # definition of computing flags
        self.CalculateReactionsFlag = True
        # Xref(n+1) = Xref(n)+dX
        self.MoveMeshFlag = True

    #
    def Initialize(self, restart):

        # definition of time scheme
        self.damp_factor_f = 0.00;
        self.damp_factor_m = -0.01;

        self.mechanical_scheme = ResidualBasedStaticScheme()
        # self.time_scheme = ResidualBasedGeneralizedAlphaScheme(self.damp_factor_f,self.damp_factor_m)
        # self.mechanical_scheme = ResidualBasedBossakScheme(self.damp_factor_m,self.dynamicstype)
        if(self.problemtype == "Mechanical"):
            self.mechanical_scheme.Check(self.model_part)

        # creating the solution strategy
        self.ReformDofSetAtEachStep = False

        if(self.ComputeMechanicsFlag):
            self.ReformDofSetAtEachStep = True

        # import solid_mechanics_python_strategy

        # self.mechanical_solver = solid_mechanics_python_strategy.ResidualStrategy(self.model_part,self.mechanical_scheme,self.solid_linear_solver,self.mechanical_convergence_criteria,self.CalculateReactionsFlag,self.ReformDofSetAtEachStep,self.LineSearch)

        self.mechanical_solver = ResidualBasedNewtonRaphsonStrategy(self.model_part, self.mechanical_scheme, self.solid_linear_solver, self.mechanical_convergence_criteria, self.max_iter, self.CalculateReactionsFlag, self.ReformDofSetAtEachStep, self.MoveMeshFlag)

        if(self.problemtype == "Mechanical"):
            self.mechanical_solver.Check();

        # problems with Reactions calculation on thermal-mechanical builder_and_solver
        # self.CalculateReactionsFlag = False
        if(restart):
            self.mechanical_solver.SetInitialized()

    #
    def Solve(self):

        if(self.ComputeMechanicsFlag):
            print(" MECHANICAL SOLUTION START ")
            (self.mechanical_solver).Solve()
            print(" MECHANICAL SOLUTION PERFORMED ")

        # move the mesh as needed
        # if(self.MoveMeshFlag == True):
        #    self.mechanical_scheme.MoveMesh(self.model_part.Nodes);
        # MoveMesh is not a function of the scheme is a function of add strategies to python
    #
    def SetEchoLevel(self, level):
        (self.mechanical_solver).SetEchoLevel(level)

    #
    def SetMaxIters(self, iters):
        self.max_iter = iters
        #(self.mechanical_solver).SetMaxIters(iters)

    #
    def Check(self):
        self.builder_and_solver.Check(self.model_part)
        self.mechanical_scheme.Check(self.model_part)
        self.mechanical_convergence_criteria.Check(self.model_part)
