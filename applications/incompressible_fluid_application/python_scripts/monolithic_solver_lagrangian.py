from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# importing the Kratos Library
from KratosMultiphysics import *
from KratosMultiphysics.IncompressibleFluidApplication import *
from KratosMultiphysics.PFEMApplication import *
from KratosMultiphysics.MeshingApplication import *
from KratosMultiphysics.ExternalSolversApplication import *

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
    model_part.AddNodalSolutionStepVariable(ADVPROJ)
    model_part.AddNodalSolutionStepVariable(DIVPROJ)
    model_part.AddNodalSolutionStepVariable(THAWONE)
    model_part.AddNodalSolutionStepVariable(THAWTWO)
    model_part.AddNodalSolutionStepVariable(REACTION)
    model_part.AddNodalSolutionStepVariable(REACTION_WATER_PRESSURE)
    model_part.AddNodalSolutionStepVariable(EXTERNAL_PRESSURE)
    model_part.AddNodalSolutionStepVariable(DISTANCE)
    model_part.AddNodalSolutionStepVariable(ARRHENIUS)
    model_part.AddNodalSolutionStepVariable(IS_WATER)

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
# self.linear_solver =  SkylineLUFactorizationSolver()
# self.linear_solver =SuperLUSolver()

        pPrecond = DiagonalPreconditioner()
# pPrecond = ILU0Preconditioner()
        self.linear_solver = BICGSTABSolver(1e-6, 5000, pPrecond)
# self.linear_solver = CGSolver(1e-6, 5000,pPrecond)

        # definition of the convergence criteria
# The argument order: VelRatioTolerance;  VelAbsTolerance;
# PrsRatioTolerance; PrsAbsTolerance;

        self.conv_criteria = VelPrCriteria(1e-7, 1e-9, 1e-7, 1e-9)

        self.dynamic_tau = 0.0
        self.oss_swith = 0
        self.max_iter = 10

        # default settings
        self.echo_level = 1
        self.CalculateReactionFlag = False
        self.ReformDofSetAtEachStep = True
        self.CalculateNormDxFlag = True
        self.MoveMeshFlag = True
        self.remeshing_flag = True

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
            self.neigh_finder = FindNodalNeighboursProcess(model_part, 20, 30)

        self.alpha_shape = 1.2
        self.h_factor = 0.4

        # assign IS_FLUID to all nodes
# for node in self.model_part.Nodes:
# node.SetSolutionStepValue(IS_FLUID,0,1.0)

        # detecting free_surface to all nodes
        for node in self.model_part.Nodes:
            if (node.GetSolutionStepValue(IS_BOUNDARY) == 1 and node.GetSolutionStepValue(IS_STRUCTURE) != 1):
                node.SetSolutionStepValue(IS_FREE_SURFACE, 0, 1.0)

        # U NEED IT FOR ALPHA-shape
        (self.neigh_finder).Execute()
        Hfinder = FindNodalHProcess(model_part)
        Hfinder.Execute()

        # runtime box
        self.box_corner1 = Array3()
        self.box_corner1[0] = box_corner1[0]
        self.box_corner1[1] = box_corner1[1]
        self.box_corner1[2] = box_corner1[2]

        self.box_corner2 = Array3()
        self.box_corner2[0] = box_corner2[0]
        self.box_corner2[1] = box_corner2[1]
        self.box_corner2[2] = box_corner2[2]

    #
    def Initialize(self, output_time_increment):
        # creating the solution strategy

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
        self.model_part.ProcessInfo.SetValue(OSS_SWITCH, self.oss_swith)

        # time increment for output
        self.output_time_increment = output_time_increment
        self.next_output_time = self.output_time_increment

#        (self.neigh_finder).Execute();

    #
    def Solve(self, time, gid_io):
# (self.neigh_finder).Execute();
# (self.solver).Solve()
# (self.solver).Clear()
# (self.PfemUtils).MarkOuterNodes(self.box_corner1,self.box_corner2,(self.model_part).Nodes );
# (self.PfemUtils).MarkExcessivelyCloseNodes((self.model_part).Nodes, .05)
# (self.node_erase_process).Execute();
# self.Remesh()
# self.OutputStep(time,gid_io)
        print("143")

        self.Remesh()
        print("145")
        (self.solver).Solve()
        # print
        # "KKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKK"
        (self.PfemUtils).MoveLonelyNodes(self.model_part)
        print("a47")
        (self.solver).Clear()
        print("149")
        self.OutputStep(time, gid_io)

    #
    def EstimateDeltaTime(self, min_dt, max_dt):
        print("Estimating delta time")
        return (
            (self.PfemUtils).EstimateDeltaTime(min_dt, max_dt, self.model_part)
        )

#    def EstimateDeltaTime(self,min_dt,max_dt):
#        print "Estimating delta time"
#        return (self.UlfUtils).EstimateDeltaTime(max_dt,domain_size)

    #
    def SetEchoLevel(self, level):
        (self.solver).SetEchoLevel(level)

#
# def Remesh(self):
#
# if (self.remeshing_flag==True):
# (self.Mesher).ReGenerateMesh("ASGS2D", "Condition2D",self.model_part,self.node_erase_process,True, True, self.alpha_shape, self.h_factor)
# (self.Mesher).ReGenerateMesh("ASGS2D", "Condition2D",self.model_part,self.node_erase_process,True, False, self.alpha_shape, self.h_factor)
#
# calculating fluid neighbours before applying boundary conditions
# (self.neigh_finder).Execute();

    #
    def Remesh(self):

        # if (self.remeshing_flag==True):
           # (self.PfemUtils).MoveLonelyNodes(self.model_part)
        (self.MeshMover).Execute()

        (self.PfemUtils).MarkOuterNodes(self.box_corner1,
                                        self.box_corner2, (self.model_part).Nodes)
        (self.PfemUtils).MarkNodesTouchingWall(
            self.model_part, self.domain_size, 0.05)
        (self.node_erase_process).Execute()

        (self.neigh_finder).ClearNeighbours()

        ((self.model_part).Elements).clear()
        ((self.model_part).Conditions).clear()

        # remesh
        if(self.domain_size == 2):
            (self.Mesher).ReGenerateMesh("ASGS2D", "Condition2D", self.model_part,
                                         self.node_erase_process, True, True, self.alpha_shape, self.h_factor)
        elif(self.domain_size == 3):
            (self.Mesher).ReGenerateMesh("ASGS3D", "Condition3D", self.model_part,
                                         self.node_erase_process, True, True, self.alpha_shape, self.h_factor)

         # calculating fluid neighbours before applying boundary conditions
        (self.neigh_finder).Execute()

        (self.PfemUtils).ApplyBoundaryConditions(self.model_part, 2)
        (self.PfemUtils).IdentifyFluidNodes(self.model_part)
        (self.PfemUtils).ApplyMinimalPressureConditions(self.model_part)

# for node in self.model_part.Nodes:
# node.SetSolutionStepValue(IS_FREE_SURFACE,0,0.0)
#
# for node in self.model_part.Nodes:
# if (node.GetSolutionStepValue(IS_BOUNDARY)==1 and node.GetSolutionStepValue(IS_STRUCTURE)!=1):
# node.SetSolutionStepValue(IS_FREE_SURFACE,0,1.0)
        #
    def FindNeighbours(self):
        (self.neigh_finder).Execute()

    #
    def OutputStep(self, time, gid_io):
        if(time >= self.next_output_time):
            self.next_output_time = self.next_output_time + \
                self.output_time_increment

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
                IS_FLUID,
                (self.model_part).Nodes,
                time,
                0)

            gid_io.Flush()
            gid_io.FinalizeResults()
    #

    def ReadRestartFile(self, FileName, nodes):
        NODES = nodes
        aaa = open(FileName)

        for line in aaa:
            print(line)
            exec(line)

    #
    def WriteRestartFile(self, FileName):
        backupfile = open(FileName + ".py", 'w')

        import restart_utilities
        restart_utilities.PrintRestart_ScalarVariable_PyFormat(
            VELOCITY_X,
            "VELOCITY_X",
            self.model_part.Nodes,
            backupfile)
        restart_utilities.PrintRestart_ScalarVariable_PyFormat(
            VELOCITY_Y,
            "VELOCITY_Y",
            self.model_part.Nodes,
            backupfile)
        restart_utilities.PrintRestart_ScalarVariable_PyFormat(
            VELOCITY_Z,
            "VELOCITY_Z",
            self.model_part.Nodes,
            backupfile)
        restart_utilities.PrintRestart_ScalarVariable_PyFormat(
            PRESSURE, "PRESSURE", self.model_part.Nodes, backupfile)
        restart_utilities.PrintRestart_ScalarVariable_PyFormat(
            DENSITY, "DENSITY", self.model_part.Nodes, backupfile)
        restart_utilities.PrintRestart_ScalarVariable_PyFormat(
            VISCOSITY, "VISCOSITY", self.model_part.Nodes, backupfile)
        restart_utilities.PrintRestart_ScalarVariable_PyFormat(
            NODAL_H, "NODAL_H", self.model_part.Nodes, backupfile)
        restart_utilities.PrintRestart_ScalarVariable_PyFormat(
            IS_STRUCTURE,
            "IS_STRUCTURE",
            self.model_part.Nodes,
            backupfile)
        restart_utilities.PrintRestart_ScalarVariable_PyFormat(
            IS_BOUNDARY,
            "IS_BOUNDARY",
            self.model_part.Nodes,
            backupfile)
        restart_utilities.PrintRestart_ScalarVariable_PyFormat(
            IS_WATER, "IS_WATER", self.model_part.Nodes, backupfile)
        restart_utilities.PrintRestart_ScalarVariable_PyFormat(
            IS_INTERFACE,
            "IS_INTERFACE",
            self.model_part.Nodes,
            backupfile)
        restart_utilities.PrintRestart_ScalarVariable_PyFormat(
            DISTANCE, "DISTANCE", self.model_part.Nodes, backupfile)
        restart_utilities.PrintRestart_ScalarVariable_PyFormat(
            DISPLACEMENT_X,
            "DISPLACEMENT_X",
            self.model_part.Nodes,
            backupfile)
        restart_utilities.PrintRestart_ScalarVariable_PyFormat(
            DISPLACEMENT_Y,
            "DISPLACEMENT_Y",
            self.model_part.Nodes,
            backupfile)
        restart_utilities.PrintRestart_ScalarVariable_PyFormat(
            DISPLACEMENT_Z,
            "DISPLACEMENT_Z",
            self.model_part.Nodes,
            backupfile)
        restart_utilities.PrintRestart_ScalarVariable_PyFormat(
            DISTANCE, "DISTANCE", self.model_part.Nodes, backupfile)
        restart_utilities.PrintRestart_ScalarVariable_PyFormat(
            ACCELERATION_X,
            "ACCELERATION_X",
            self.model_part.Nodes,
            backupfile)
        restart_utilities.PrintRestart_ScalarVariable_PyFormat(
            ACCELERATION_Y,
            "ACCELERATION_Y",
            self.model_part.Nodes,
            backupfile)
        restart_utilities.PrintRestart_ScalarVariable_PyFormat(
            ACCELERATION_Z,
            "ACCELERATION_Z",
            self.model_part.Nodes,
            backupfile)
        restart_utilities.PrintRestart_ScalarVariable_PyFormat(
            MESH_VELOCITY_X,
            "MESH_VELOCITY_X",
            self.model_part.Nodes,
            backupfile)
        restart_utilities.PrintRestart_ScalarVariable_PyFormat(
            MESH_VELOCITY_Y,
            "MESH_VELOCITY_Y",
            self.model_part.Nodes,
            backupfile)
        restart_utilities.PrintRestart_ScalarVariable_PyFormat(
            MESH_VELOCITY_Z,
            "MESH_VELOCITY_Z",
            self.model_part.Nodes,
            backupfile)
        restart_utilities.PrintRestart_ScalarVariable_PyFormat(
            BODY_FORCE_X,
            "BODY_FORCE_X",
            self.model_part.Nodes,
            backupfile)
        restart_utilities.PrintRestart_ScalarVariable_PyFormat(
            BODY_FORCE_Y,
            "BODY_FORCE_Y",
            self.model_part.Nodes,
            backupfile)
        restart_utilities.PrintRestart_ScalarVariable_PyFormat(
            BODY_FORCE_Z,
            "BODY_FORCE_Z",
            self.model_part.Nodes,
            backupfile)

        restart_utilities.PrintRestartFixity_PyFormat(
            VELOCITY_X,
            "VELOCITY_X",
            self.model_part.Nodes,
            backupfile)
        restart_utilities.PrintRestartFixity_PyFormat(
            VELOCITY_Y,
            "VELOCITY_Y",
            self.model_part.Nodes,
            backupfile)
        restart_utilities.PrintRestartFixity_PyFormat(
            VELOCITY_Z,
            "VELOCITY_Z",
            self.model_part.Nodes,
            backupfile)
        restart_utilities.PrintRestartFixity_PyFormat(
            PRESSURE,
            "PRESSURE",
            self.model_part.Nodes,
            backupfile)

        # restart_utilities.PrintRestart_Position_PyFormat(self.model_part.Nodes,backupfile)

        backupfile.close()
    #
