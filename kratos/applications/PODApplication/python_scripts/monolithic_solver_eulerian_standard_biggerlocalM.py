from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# -*- coding: utf-8 -*-
# importing the Kratos Library
from Kratos import *
from KratosIncompressibleFluidApplication import *
from KratosPodApplication import *
# from KratosExternalSolversApplication import *
# from KratosStructuralApplication import *


def AddVariables(model_part):
    model_part.AddNodalSolutionStepVariable(VELOCITY);
    model_part.AddNodalSolutionStepVariable(ACCELERATION);
    model_part.AddNodalSolutionStepVariable(MESH_VELOCITY);
    model_part.AddNodalSolutionStepVariable(PRESSURE);
    model_part.AddNodalSolutionStepVariable(AIR_PRESSURE);
    model_part.AddNodalSolutionStepVariable(IS_FLUID);
    model_part.AddNodalSolutionStepVariable(IS_POROUS);
    model_part.AddNodalSolutionStepVariable(IS_STRUCTURE);
    model_part.AddNodalSolutionStepVariable(IS_FREE_SURFACE);
    model_part.AddNodalSolutionStepVariable(IS_INTERFACE);
    model_part.AddNodalSolutionStepVariable(IS_BOUNDARY);
    model_part.AddNodalSolutionStepVariable(DISPLACEMENT);
    model_part.AddNodalSolutionStepVariable(VISCOSITY);
    model_part.AddNodalSolutionStepVariable(DENSITY);
    model_part.AddNodalSolutionStepVariable(POROSITY);
    model_part.AddNodalSolutionStepVariable(DENSITY_AIR);
    model_part.AddNodalSolutionStepVariable(AIR_SOUND_VELOCITY);
    model_part.AddNodalSolutionStepVariable(SOUND_VELOCITY);
    model_part.AddNodalSolutionStepVariable(BODY_FORCE);
    model_part.AddNodalSolutionStepVariable(NODAL_AREA);
    model_part.AddNodalSolutionStepVariable(NODAL_H);
    model_part.AddNodalSolutionStepVariable(ADVPROJ);
    model_part.AddNodalSolutionStepVariable(DIVPROJ);
    model_part.AddNodalSolutionStepVariable(THAWONE);
    model_part.AddNodalSolutionStepVariable(THAWTWO);
    model_part.AddNodalSolutionStepVariable(REACTION);
    model_part.AddNodalSolutionStepVariable(REACTION_WATER_PRESSURE);
    model_part.AddNodalSolutionStepVariable(EXTERNAL_PRESSURE);
    model_part.AddNodalSolutionStepVariable(WATER_PRESSURE);
    model_part.AddNodalSolutionStepVariable(AIR_PRESSURE_DT);
    model_part.AddNodalSolutionStepVariable(ARRHENIUS);
    model_part.AddNodalSolutionStepVariable(FLAG_VARIABLE);
    model_part.AddNodalSolutionStepVariable(NORMAL);

    print("variables for the dynamic structural solution added correctly")


def AddDofs(model_part):
    for node in model_part.Nodes:
        # adding dofs
        node.AddDof(VELOCITY_X, REACTION_X);
        node.AddDof(VELOCITY_Y, REACTION_Y);
        node.AddDof(VELOCITY_Z, REACTION_Z);
        node.AddDof(PRESSURE, REACTION_WATER_PRESSURE);
        node.AddDof(AIR_PRESSURE, REACTION_AIR_PRESSURE);

    print("dofs for the monolithic solver added correctly")


class MonolithicSolver:
    #

    def __init__(self, model_part, domain_size):

        self.model_part = model_part

        self.alpha = -0.3
        self.move_mesh_strategy = 0
        self.time_scheme = ResidualBasedPredictorCorrectorVelocityBossakScheme(self.alpha, self.move_mesh_strategy)
        # definition of the solvers
        self.linear_solver = SkylineLUFactorizationSolver()
# self.linear_solver =SuperLUSolver()

        # pPrecond = DiagonalPreconditioner()
# pPrecond = ILU0Preconditioner()
        # self.linear_solver =  BICGSTABSolver(1e-6, 5000,pPrecond)

        # definition of the convergence criteria
# The argument order: VelRatioTolerance;	VelAbsTolerance; PrsRatioTolerance; PrsAbsTolerance;
        self.conv_criteria = UPCriteria(1e-7, 1e-7, 1e-3, 1e-7)
       # self.conv_criteria = UPCriteria(1e-12,1e-14,1e-15,1e-17)
        self.model_part.ProcessInfo.SetValue(DYNAMIC_TAU, 0.001);

        self.max_iter = 100

        # default settings
        self.echo_level = 0
        self.CalculateReactionFlag = True
        self.ReformDofSetAtEachStep = True
        self.CalculateNormDxFlag = True
        self.MoveMeshFlag = False

# print "Construction monolithic solver finished"

    #
    def Initialize(self):
        self.builder_and_solver = ResidualBasedEliminationBuilderAndSolverStandard_BiggerLocalM(self.linear_solver)

        # Note that the strategy asks for a solver but doesn't use it (when
        # called using this constructor). This is good, as this builder and
        # solver uses 2 different solvers and one of them is given here
        # arbitrarily
        self.solver = ResidualBasedNewtonRaphsonStrategy\
            (self.model_part,
            self.time_scheme,
            self.linear_solver,
            self.conv_criteria,
            self.builder_and_solver,
            self.max_iter,
            self.CalculateReactionFlag,
            self.ReformDofSetAtEachStep,
                       self.MoveMeshFlag)

        (self.solver).SetEchoLevel(self.echo_level)

        # creating the solution strategy

        # self.solver = ResidualBasedNewtonRaphsonStrategy(self.model_part,self.time_scheme,self.linear_solver,self.conv_criteria,self.max_iter,self.CalculateReactionFlag, self.ReformDofSetAtEachStep,self.MoveMeshFlag)
        #(self.solver).SetEchoLevel(self.echo_level)
# print "Initialization monolithic solver finished"

    #
    def Solve(self):
# print "*****************entering solve?????????????"
        (self.solver).Solve()
# print "solving step monolithic solver finished"

    #
    def SetEchoLevel(self, level):
        (self.solver).SetEchoLevel(level)

    #
