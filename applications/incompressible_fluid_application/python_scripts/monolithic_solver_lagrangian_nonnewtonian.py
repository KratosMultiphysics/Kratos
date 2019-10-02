from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# importing the Kratos Library
from KratosMultiphysics import *
from KratosMultiphysics.IncompressibleFluidApplication import *
from KratosMultiphysics.PFEMApplication import *
from KratosMultiphysics.MeshingApplication import *
from KratosMultiphysics.ExternalSolversApplication import *

import math


def AddVariables(model_part):
    model_part.AddNodalSolutionStepVariable(VELOCITY)
    model_part.AddNodalSolutionStepVariable(ACCELERATION)
    model_part.AddNodalSolutionStepVariable(MESH_VELOCITY)
    model_part.AddNodalSolutionStepVariable(PRESSURE)
    model_part.AddNodalSolutionStepVariable(AIR_PRESSURE)
    model_part.AddNodalSolutionStepVariable(WATER_PRESSURE)
    model_part.AddNodalSolutionStepVariable(IS_FLUID)
    model_part.AddNodalSolutionStepVariable(IS_POROUS)
    model_part.AddNodalSolutionStepVariable(IS_STRUCTURE)
    model_part.AddNodalSolutionStepVariable(IS_FREE_SURFACE)
    model_part.AddNodalSolutionStepVariable(IS_INTERFACE)
    model_part.AddNodalSolutionStepVariable(IS_BOUNDARY)
    model_part.AddNodalSolutionStepVariable(DISPLACEMENT)
    model_part.AddNodalSolutionStepVariable(VISCOSITY)
    model_part.AddNodalSolutionStepVariable(DENSITY)
    model_part.AddNodalSolutionStepVariable(DENSITY_AIR)
    model_part.AddNodalSolutionStepVariable(AIR_SOUND_VELOCITY)
    model_part.AddNodalSolutionStepVariable(SOUND_VELOCITY)
    model_part.AddNodalSolutionStepVariable(BODY_FORCE)
    model_part.AddNodalSolutionStepVariable(NODAL_AREA)
    model_part.AddNodalSolutionStepVariable(NODAL_H)
    model_part.AddNodalSolutionStepVariable(NORMAL)
    model_part.AddNodalSolutionStepVariable(ADVPROJ)
    model_part.AddNodalSolutionStepVariable(DIVPROJ)
    model_part.AddNodalSolutionStepVariable(NORMAL)
    model_part.AddNodalSolutionStepVariable(THAWONE)
    model_part.AddNodalSolutionStepVariable(THAWTWO)
    model_part.AddNodalSolutionStepVariable(REACTION)
    model_part.AddNodalSolutionStepVariable(REACTION_WATER_PRESSURE)
    model_part.AddNodalSolutionStepVariable(EXTERNAL_PRESSURE)
    model_part.AddNodalSolutionStepVariable(DISTANCE)
    model_part.AddNodalSolutionStepVariable(ARRHENIUS)
    model_part.AddNodalSolutionStepVariable(IS_WATER)
    model_part.AddNodalSolutionStepVariable(INTERNAL_FRICTION_ANGLE)
    model_part.AddNodalSolutionStepVariable(TAU)
    model_part.AddNodalSolutionStepVariable(MU)
    model_part.AddNodalSolutionStepVariable(YIELD_STRESS)
    model_part.AddNodalSolutionStepVariable(EQ_STRAIN_RATE)
    # WATER PRESSURE GRADIENT
    model_part.AddNodalSolutionStepVariable(PRESS_PROJ)

    print("variables for the dynamic structural solution added correctly")


def AddDofs(model_part):
    for node in model_part.Nodes:
        # adding dofs
        node.AddDof(VELOCITY_X, REACTION_X)
        node.AddDof(VELOCITY_Y, REACTION_Y)
        node.AddDof(VELOCITY_Z, REACTION_Z)
        node.AddDof(PRESSURE, REACTION_WATER_PRESSURE)
        node.AddDof(AIR_PRESSURE, REACTION_AIR_PRESSURE)

    print("dofs for the monolithic solver added correctly")


class MonolithicSolver:
    #

    def __init__(self, model_part, domain_size, box_corner1, box_corner2):

        self.model_part = model_part
        self.domain_size = domain_size
        self.alpha = -0.1
        self.move_mesh_strategy = 2
        self.time_scheme = ResidualBasedPredictorCorrectorVelocityBossakScheme(
            self.alpha, self.move_mesh_strategy)
        # definition of the solvers
# self.linear_solver = SkylineLUFactorizationSolver()
        self.linear_solver = SuperLUSolver()
# self.linear_solver = SuperLUIterativeSolver()


# pPrecond = DiagonalPreconditioner()
# pPrecond = ILU0Preconditioner()
# self.linear_solver =  BICGSTABSolver(1e-6, 5000,pPrecond)
# self.linear_solver = CGSolver(1e-6, 5000,pPrecond)

        # definition of the convergence criteria
# The argument order: VelRatioTolerance;  VelAbsTolerance;
# PrsRatioTolerance; PrsAbsTolerance;
        self.conv_criteria = VelPrCriteria(5e-4, 1e-5, 5e-3, 1e-5)

        self.max_iter = 20

        # default settings
        self.dynamic_tau = 0.0
        self.oss_switch = 0
        # m_coeff in exponencial viscosity model
        self.regularization_coef = 3000
        self.echo_level = 1
        self.CalculateReactionFlag = False
        self.ReformDofSetAtEachStep = True
        self.CalculateNormDxFlag = True
        self.MoveMeshFlag = True
# self.remeshing_flag = True
        # MESH CHANGES
      #  self.mark_close_nodes_process = MarkCloseNodesProcess(model_part);
        self.PfemUtils = PfemUtils()
        self.MeshMover = MoveMeshProcess(self.model_part)
        self.node_erase_process = NodeEraseProcess(model_part)
        if(domain_size == 2):
            self.Mesher = TriGenPFEMModeler()
            self.neigh_finder = FindNodalNeighboursProcess(model_part, 9, 18)
        elif(domain_size == 3):
            self.Mesher = TetGenPfemModeler()
# self.Mesher = TetGenPfemRefineFace()
            self.neigh_finder = FindNodalNeighboursProcess(model_part, 20, 30)

        self.alpha_shape = 1.6
        self.h_factor = 0.4
# self.h_factor = 0.7
        # detecting free_surface to all nodes
        for node in self.model_part.Nodes:
            if (node.GetSolutionStepValue(IS_BOUNDARY) == 1 and node.GetSolutionStepValue(IS_STRUCTURE) != 1):
                node.SetSolutionStepValue(IS_FREE_SURFACE, 0, 1.0)

        # U NEED IT FOR ALPHA-shape
        (self.neigh_finder).Execute()
        Hfinder = FindNodalHProcess(model_part)
        Hfinder.Execute()

        # runtime box
        self.box_corner1 = box_corner1
        self.box_corner2 = box_corner2

    #
    def Initialize(self, output_time_increment):
        # creating the solution strategy
        print(self.linear_solver)

        self.solver = ResidualBasedNewtonRaphsonStrategy(
            self.model_part,
            self.time_scheme,
            self.linear_solver,
            self.conv_criteria,
            self.max_iter,
            self.CalculateReactionFlag,
            self.ReformDofSetAtEachStep,
            self.MoveMeshFlag)
        (self.solver).SetEchoLevel(self.echo_level)

        self.model_part.ProcessInfo.SetValue(DYNAMIC_TAU, self.dynamic_tau)
        self.model_part.ProcessInfo.SetValue(OSS_SWITCH, self.oss_switch)
        self.model_part.ProcessInfo.SetValue(M, self.regularization_coef)
# self.model_part.ProcessInfo.SetValue(SCALE, self.regularization_coef );

        # time increment for output
        self.output_time_increment = output_time_increment
        self.next_output_time = self.output_time_increment

#        (self.neigh_finder).Execute();

        print("Remesh performed passing the element and the condition directly - not their name")
        self.reference_element = self.model_part.Elements[1]
# print "refenence element:   ",self.reference_element
        self.reference_condition = self.model_part.Conditions[1]
# print "refenence condition:   ",self.reference_condition
    #

    def Solve(self, time, gid_io):

        self.Remesh()

        (self.solver).Solve()

# self.RestoreOldPosition()
        (self.PfemUtils).MoveLonelyNodes(self.model_part)

        (self.solver).Clear()

        self.OutputStep(time, gid_io)

    #
    def EstimateDeltaTime(self, min_dt, max_dt):
        print("Estimating delta time")
        return (
            (self.PfemUtils).EstimateDeltaTime(min_dt, max_dt, self.model_part)
        )

    #
    def SetEchoLevel(self, level):
        (self.solver).SetEchoLevel(level)

    #
    def Remesh(self):

       # (self.PfemUtils).MoveLonelyNodes(self.model_part)
        (self.MeshMover).Execute()

        (self.PfemUtils).MarkOuterNodes(self.box_corner1,
                                        self.box_corner2, (self.model_part).Nodes)
        if(self.domain_size == 2):
            (self.PfemUtils).MarkNodesTouchingWall(
                self.model_part, self.domain_size, 0.08)
        if(self.domain_size == 3):
            (self.PfemUtils).MarkNodesTouchingWall(
                self.model_part, self.domain_size, 0.1)
# (self.ActOnWalls).Execute();
        (self.node_erase_process).Execute()

# if (self.remeshing_flag==True):
        (self.neigh_finder).ClearNeighbours()

        ((self.model_part).Elements).clear()
        ((self.model_part).Conditions).clear()

        # remesh
        if(self.domain_size == 2):
# (self.Mesher).ReGenerateMesh("NoNewtonianASGS2D", "Condition2D",self.model_part,self.node_erase_process,True, True, self.alpha_shape, self.h_factor)
# (self.Mesher).ReGenerateMesh("BinghamNonNewtonianASGS2D", "Condition2D",self.model_part,self.node_erase_process,True, True, self.alpha_shape, self.h_factor)
            (self.Mesher).ReGenerateMesh(self.model_part, self.reference_element,
                                         self.reference_condition, self.node_erase_process, True, True, self.alpha_shape, self.h_factor)
        elif(self.domain_size == 3):
# (self.Mesher).ReGenerateMesh("NoNewtonianASGS3D", "Condition3D",self.model_part,self.node_erase_process,True, True, self.alpha_shape, self.h_factor)
# (self.Mesher).ReGenerateMesh("BinghamNonNewtonianASGS2D", "Condition3D",self.model_part,self.node_erase_process,True, True, self.alpha_shape, self.h_factor)
            (self.Mesher).ReGenerateMesh(self.model_part, self.reference_element,
                                         self.reference_condition, self.node_erase_process, True, True, self.alpha_shape, self.h_factor)

        print("regenerated mesh")

         # calculating fluid neighbours before applying boundary conditions
        (self.neigh_finder).Execute()
        print("found neighbours")
        (self.PfemUtils).ApplyBoundaryConditions(self.model_part, 2)
        (self.PfemUtils).IdentifyFluidNodes(self.model_part)
        (self.PfemUtils).ApplyMinimalPressureConditions(self.model_part)
        print("applied BC")

# for node in self.model_part.Nodes:
# node.SetSolutionStepValue(IS_FREE_SURFACE,0,0.0)
#
# for node in self.model_part.Nodes:
# if (node.GetSolutionStepValue(IS_BOUNDARY)==1 and node.GetSolutionStepValue(IS_STRUCTURE)!=1):
# node.SetSolutionStepValue(IS_FREE_SURFACE,0,1.0)

    #
# def Remesh3D(self, param):
# self.remeshing_flag==False
# delta_displ_max = 0.0
# nodal_h_max = 0.0
# for node in self.model_part.Nodes:
# displacement is the total displacement from the beginning of the simulation.
# We have to consider the Delta_displacement
# displX = node.GetSolutionStepValue(DISPLACEMENT_X)
# displY = node.GetSolutionStepValue(DISPLACEMENT_Y)
# displZ = node.GetSolutionStepValue(DISPLACEMENT_Z)
# old_displX = node.GetSolutionStepValue(DISPLACEMENT_X,1)
# old_displY = node.GetSolutionStepValue(DISPLACEMENT_Y,1)
# old_displZ = node.GetSolutionStepValue(DISPLACEMENT_Z,1)
#
# displ = math.sqrt(displX*displX + displY*displY + displZ*displZ)
# delta_displ_square = (displX - old_displX)*(displX - old_displX) + (displY - old_displY)*(displY - old_displY) + (displZ - old_displZ)*(displZ - old_displZ)
# if (delta_displ_square > delta_displ_max):#value independent from discretization----> not good to compare with mesh dimension....
# delta_displ_max = delta_displ_square
# delta_displ_max = math.sqrt(delta_displ_max)
# if(node.GetSolutionStepValue(NODAL_H) > nodal_h_max):#aprox max mesh dimension
# nodal_h_max = node.GetSolutionStepValue(NODAL_H)
#
#
# if(delta_displ_max > param):
# self.remeshing_flag==True
#
# (self.MeshMover).Execute();
#
# (self.PfemUtils).MarkOuterNodes(self.box_corner1,self.box_corner2,(self.model_part).Nodes );
# (self.PfemUtils).MarkNodesTouchingWall(self.model_part, self.domain_size, 0.08)
#
# (self.node_erase_process).Execute();
#
# if(self.remeshing_flag==True):
# print Time, "-------------------------------------------------REMESH_3D-----------"
#
# (self.neigh_finder).ClearNeighbours();
#
# ((self.model_part).Elements).clear();
# ((self.model_part).Conditions).clear();
#
# remesh
# if(self.domain_size == 2):
# (self.Mesher).ReGenerateMesh("NoNewtonianASGS2D", "Condition2D",self.model_part,self.node_erase_process,True, True, self.alpha_shape, self.h_factor)
# elif(self.domain_size == 3):
# (self.Mesher).ReGenerateMesh("NoNewtonianASGS3D", "Condition3D",self.model_part,self.node_erase_process,True, True, self.alpha_shape, self.h_factor)
#
# print "regenerated mesh"
#
# calculating fluid neighbours before applying boundary conditions
# (self.neigh_finder).Execute();
# print "found neighbours"
# (self.PfemUtils).ApplyBoundaryConditions(self.model_part,2);
# (self.PfemUtils).IdentifyFluidNodes(self.model_part);
# (self.PfemUtils).ApplyMinimalPressureConditions(self.model_part);
# print "applied BC"
#
# self.remeshing_flag==False

    #
    def FindNeighbours(self):
        (self.neigh_finder).Execute()

    #
    def RestoreOldPosition(self):

        for node in self.model_part.Nodes:
            # displacement is the total displacement from the beginning of the simulation.
            # We have to consider the Delta_displacement
            displX = node.GetSolutionStepValue(DISPLACEMENT_X)
            displY = node.GetSolutionStepValue(DISPLACEMENT_Y)
            displZ = node.GetSolutionStepValue(DISPLACEMENT_Z)
            old_displX = node.GetSolutionStepValue(DISPLACEMENT_X, 1)
            old_displY = node.GetSolutionStepValue(DISPLACEMENT_Y, 1)
            old_displZ = node.GetSolutionStepValue(DISPLACEMENT_Z, 1)

            displ = math.sqrt(
                displX *
                displX +
                displY *
                displY +
                displZ *
                displZ)
            delta_displ = (
                displX - old_displX) * (
                    displX - old_displX) + (
                        displY - old_displY) * (
                            displY - old_displY) + (
                                displZ - old_displZ) * (
                                    displZ - old_displZ)
            delta_displ = math.sqrt(delta_displ)

            # 1cm for the moment but this parameter should be passed by the
            # user...
            if(delta_displ < 0.001):
                node.SetSolutionStepValue(DISPLACEMENT_X, 0, old_displX)
                node.SetSolutionStepValue(DISPLACEMENT_Y, 0, old_displY)
                node.SetSolutionStepValue(DISPLACEMENT_Z, 0, old_displZ)

        # to update the position with the new displacement
        (self.MeshMover).Execute()

    #
    def OutputStep(self, time, gid_io):
        if(time >= self.next_output_time):
            self.next_output_time = self.next_output_time + \
                self.output_time_increment
# print "output time increment", self.output_time_increment

            # writing mesh
            gid_io.InitializeMesh(time)
            gid_io.WriteNodeMesh((self.model_part).GetMesh())
            gid_io.WriteMesh((self.model_part).GetMesh())
            gid_io.FinalizeMesh()

            gid_io.InitializeResults(time, (self.model_part).GetMesh())

            gid_io.WriteNodalResults(
                PRESSURE,
                (self.model_part).Nodes,
                time,
                0)
            gid_io.WriteNodalResults(
                IS_FREE_SURFACE,
                (self.model_part).Nodes,
                time,
                0)
            gid_io.WriteNodalResults(
                IS_BOUNDARY,
                (self.model_part).Nodes,
                time,
                0)
            gid_io.WriteNodalResults(
                IS_STRUCTURE,
                (self.model_part).Nodes,
                time,
                0)
            gid_io.WriteNodalResults(
                VELOCITY,
                (self.model_part).Nodes,
                time,
                0)
            gid_io.WriteNodalResults(
                MESH_VELOCITY,
                (self.model_part).Nodes,
                time,
                0)
            gid_io.WriteNodalResults(
                DENSITY,
                (self.model_part).Nodes,
                time,
                0)
            gid_io.WriteNodalResults(
                VISCOSITY,
                (self.model_part).Nodes,
                time,
                0)
            gid_io.WriteNodalResults(
                BODY_FORCE,
                (self.model_part).Nodes,
                time,
                0)
            gid_io.WriteNodalResults(
                DISPLACEMENT,
                (self.model_part).Nodes,
                time,
                0)
            gid_io.WriteNodalResults(
                INTERNAL_FRICTION_ANGLE,
                (self.model_part).Nodes,
                time,
                0)
            gid_io.WriteNodalResults(
                YIELD_STRESS,
                (self.model_part).Nodes,
                time,
                0)
            gid_io.PrintOnGaussPoints(EQ_STRAIN_RATE, self.model_part, time)
            gid_io.PrintOnGaussPoints(MU, self.model_part, time)
            gid_io.PrintOnGaussPoints(TAU, self.model_part, time)
# gid_io.WriteNodalResults(EFFECTIVE_VISCOSITY, (self.model_part).Nodes, time, 0);
            gid_io.WriteNodalResults(
                IS_FLUID,
                (self.model_part).Nodes,
                time,
                0)
            gid_io.WriteNodalResults(
                NODAL_H,
                (self.model_part).Nodes,
                time,
                0)
            gid_io.WriteNodalResults(
                PRESS_PROJ,
                (self.model_part).Nodes,
                time,
                0)
            # WATER PRESSURE GRADIENT

            gid_io.Flush()
            gid_io.FinalizeResults()
