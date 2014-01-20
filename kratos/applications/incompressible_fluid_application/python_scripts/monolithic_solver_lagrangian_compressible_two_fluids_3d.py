from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# importing the Kratos Library
from KratosMultiphysics import *
from KratosMultiphysics.IncompressibleFluidApplication import *
from KratosMultiphysics.PFEMApplicatio import *
from KratosMultiphysics.MeshingApplication import *
from KratosMultiphysics.ExternalSolversApplication import *
# Check that KratosMultiphysics was imported in the main script
CheckForPreviousImport()


def AddVariables(model_part):
    model_part.AddNodalSolutionStepVariable(VELOCITY)
    model_part.AddNodalSolutionStepVariable(ACCELERATION)
    model_part.AddNodalSolutionStepVariable(MESH_VELOCITY)
    model_part.AddNodalSolutionStepVariable(PRESSURE)
    model_part.AddNodalSolutionStepVariable(AIR_PRESSURE)
    model_part.AddNodalSolutionStepVariable(WATER_PRESSURE)
    model_part.AddNodalSolutionStepVariable(AIR_PRESSURE_DT)
    model_part.AddNodalSolutionStepVariable(WATER_PRESSURE_DT)
    model_part.AddNodalSolutionStepVariable(IS_FLUID)
    model_part.AddNodalSolutionStepVariable(IS_WATER)
    model_part.AddNodalSolutionStepVariable(IS_VISITED)
    model_part.AddNodalSolutionStepVariable(IS_POROUS)
    model_part.AddNodalSolutionStepVariable(IS_STRUCTURE)
    model_part.AddNodalSolutionStepVariable(IS_FREE_SURFACE)
    model_part.AddNodalSolutionStepVariable(IS_INTERFACE)
    model_part.AddNodalSolutionStepVariable(IS_BOUNDARY)
    model_part.AddNodalSolutionStepVariable(ERASE_FLAG)
    model_part.AddNodalSolutionStepVariable(DISPLACEMENT)
    model_part.AddNodalSolutionStepVariable(VISCOSITY)
    model_part.AddNodalSolutionStepVariable(VISCOSITY_AIR)
    model_part.AddNodalSolutionStepVariable(VISCOSITY_WATER)
    model_part.AddNodalSolutionStepVariable(DENSITY)
    model_part.AddNodalSolutionStepVariable(DENSITY_AIR)
    model_part.AddNodalSolutionStepVariable(DENSITY_WATER)
    model_part.AddNodalSolutionStepVariable(AIR_SOUND_VELOCITY)
    model_part.AddNodalSolutionStepVariable(WATER_SOUND_VELOCITY)
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
    model_part.AddNodalSolutionStepVariable(ARRHENIUS)
    model_part.AddNodalSolutionStepVariable(DISTANCE)
    model_part.AddNodalSolutionStepVariable(AUX_INDEX)
    model_part.AddNodalSolutionStepVariable(NORMAL)
    model_part.AddNodalSolutionStepVariable(RHS)
    model_part.AddNodalSolutionStepVariable(RHS_WATER)
    model_part.AddNodalSolutionStepVariable(RHS_AIR)
    model_part.AddNodalSolutionStepVariable(NODAL_MASS)
    model_part.AddNodalSolutionStepVariable(NODAL_PAUX)
    model_part.AddNodalSolutionStepVariable(NODAL_MAUX)

    print("variables for monolithic solver lagrangian compressible 3D solution added correctly")


def AddDofs(model_part):
    for node in model_part.Nodes:
        # adding dofs
        node.AddDof(VELOCITY_X, REACTION_X)
        node.AddDof(VELOCITY_Y, REACTION_Y)
        node.AddDof(VELOCITY_Z, REACTION_Z)
        node.AddDof(WATER_PRESSURE, REACTION_WATER_PRESSURE)
        node.AddDof(AIR_PRESSURE, REACTION_AIR_PRESSURE)

    print("variables for monolithic solver lagrangian compressible 3D solution added correctly")


class MonolithicSolver:
    #

    def __init__(self, model_part, domain_size, box_corner1, box_corner2):

        self.model_part = model_part

        self.alpha = -0.1
        self.move_mesh_strategy = 2
        self.time_scheme = ResidualBasedPredictorCorrectorVelocityBossakSchemeCompressible(
            self.alpha, self.move_mesh_strategy)
        # definition of the solvers
# self.linear_solver =  SkylineLUFactorizationSolver()
# self.linear_solver =SuperLUSolver()

        pPrecond = DiagonalPreconditioner()
# pPrecond = ILU0Preconditioner()
        self.linear_solver = BICGSTABSolver(1e-6, 5000, pPrecond)
# self.linear_solver = CGSolver(1e-6, 5000,pPrecond)

        # definition of the convergence criteria
        self.conv_criteria = UPCriteria(1e-7, 1e-9, 1e-7, 1e-9)
       # self.conv_criteria = UPCriteria(1e-12,1e-14,1e-15,1e-17)

        self.max_iter = 30

        self.SetDivided = ElemBasedBCUtilities(model_part)
        self.ChooseElement = ChooseElementProcess(
            model_part,
            3,
            "ASGSCOMPPRDC3D",
            "ASGSCompressible3D")
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

        self.node_erase_process = NodeEraseProcess(self.model_part)

        self.Mesher = TetGenPfemModeler()
       # self.Mesher = TetGenPfemRefineFace()

        self.neigh_finder = FindNodalNeighboursProcess(model_part, 9, 18)
        self.elem_neighbor_finder = FindElementalNeighboursProcess(
            model_part, 3, 20)

        self.alpha_shape = 1000000.0
        self.h_factor = 0.5
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

        # find elem neighbors for tube
        (self.elem_neighbor_finder).ClearNeighbours()
        (self.elem_neighbor_finder).Execute()

        # runtime box
        self.box_corner1 = box_corner1
        self.box_corner2 = box_corner2

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

        # time increment for output
        self.output_time_increment = output_time_increment
        self.next_output_time = self.output_time_increment

#        (self.neigh_finder).Execute();

    #
    def Solve(self, time, gid_io):
        (self.neigh_finder).Execute()
# (self.solver).Solve()
# (self.solver).Clear()
# (self.PfemUtils).MarkOuterNodes(self.box_corner1,self.box_corner2,(self.model_part).Nodes );
# (self.PfemUtils).MarkExcessivelyCloseNodes((self.model_part).Nodes, .05)
# (self.node_erase_process).Execute();
# self.Remesh()
# self.OutputStep(time,gid_io)
       # self.CalculateDistanceAndDiviedSet(3);
        print("143")
       # self.DistToH()
       # self.Remesh()
        print("145")
        (self.solver).Solve()
        print("a47")
        (self.solver).Clear()
        print("149")
        self.OutputStep(time, gid_io)

    #
    def EstimateDeltaTime(self, min_dt, max_dt):
        print("Estimating delta time")
        calc_dt = (
            self.PfemUtils).EstimateDeltaTime(
                min_dt,
                max_dt,
                self.model_part)
        # cfl_dt=(self.PfemUtils).CFLdeltaT(1.0,max_dt,self.model_part)
        # max_dt = cfl_dt
        # print"CFL_CHOICE",cfl_dt
        # calc_dt=(self.PfemUtils).ExactDtEstimate(max_dt,self.model_part)

       # print "calculated dt"
        return calc_dt

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

        if (self.remeshing_flag):
        #    self.CalculateDistanceAndDiviedSet(3);
            (self.PfemUtils).MoveLonelyNodes(self.model_part)
            #(self.MeshMover).Execute();
            print(self.box_corner1)
            (self.PfemUtils).MarkOuterNodes(
                self.box_corner1, self.box_corner2, (self.model_part).Nodes)
            (self.PfemUtils).MarkNodesTouchingWall(self.model_part, 3, .05)
            print("after nodes touching wall")
            (self.PfemUtils).MarkExcessivelyCloseNodes(
                (self.model_part).Nodes, 0.5)
            print("after excessively close nodes")
            (self.PfemUtils).MarkNodesTouchingInterface(self.model_part, 3, .1)
            print("after MarkNodesTouchingInterface")
# FIND NEIGHBOUR ELEMENTS AND COLORing
            (self.elem_neighbor_finder).ClearNeighbours()
            (self.elem_neighbor_finder).Execute()
            (self.PfemUtils).ColourAirWaterElement(self.model_part, 3)

#
# (self.PfemUtils).InterfaceDetecting(self.model_part,2, .9)
# (self.PfemUtils).ChangeWallWaterFlag(self.model_part,2)
# (self.PfemUtils).ChangeInterfaceWaterFlag(self.model_part,2)

# for node in (self.model_part).Nodes:
# if(node.GetSolutionStepValue(IS_INTERFACE) == 1.0):
# print node.GetValue(ERASE_FLAG)

            #(self.node_erase_process).Execute(); to be able to compute neighbors earase process is done inside the mesher

            (self.neigh_finder).ClearNeighbours()
            (self.neigh_finder).Execute()

            (self.Mesher).ReGenerateMesh("ASGSCompressible3D", "Condition3D", self.model_part,
                                         self.node_erase_process, True, True, self.alpha_shape, self.h_factor)
            #(self.Mesher).ReGenerateMesh("ASGSCompressible3D", "Condition3D",self.model_part,(self.structure_model_part).Elements,self.node_erase_process,True, True, self.alpha_shape, self.h_factor)

# (self.Mesher).ReGenerateMesh("ASGS2D", "Condition2D",self.model_part,self.node_erase_process,True, False, self.alpha_shape, self.h_factor)
            print("after remesh")
             # calculating fluid neighbours before applying boundary
             # conditionsTrue
            (self.elem_neighbor_finder).ClearNeighbours()
            (self.elem_neighbor_finder).Execute()
            (self.neigh_finder).ClearNeighbours()
            (self.neigh_finder).Execute()
            print("neighbors are calculated")

            (self.PfemUtils).ColourAirWaterElement(self.model_part, 3)
            (self.PfemUtils).InterfaceDetecting(self.model_part, 3, .9)
            (self.ChooseElement).Execute()

             # calculating fluid neighbours before applying boundary conditions
        #    (self.neigh_finder).ClearNeighbours();
         #   (self.neigh_finder).Execute();
            (self.PfemUtils).ApplyBoundaryConditions(self.model_part, 3)
            (self.PfemUtils).IdentifyFluidNodes(self.model_part)

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
                EXTERNAL_PRESSURE,
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
                IS_INTERFACE,
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
                AIR_PRESSURE,
                (self.model_part).Nodes,
                time,
                0)
            gid_io.WriteNodalResults(
                WATER_PRESSURE,
                (self.model_part).Nodes,
                time,
                0)
            gid_io.WriteNodalResults(
                DENSITY_AIR,
                (self.model_part).Nodes,
                time,
                0)
            gid_io.WriteNodalResults(
                DENSITY_WATER,
                (self.model_part).Nodes,
                time,
                0)
            gid_io.WriteNodalResults(
                AIR_SOUND_VELOCITY,
                (self.model_part).Nodes,
                time,
                0)
            gid_io.WriteNodalResults(
                WATER_SOUND_VELOCITY,
                (self.model_part).Nodes,
                time,
                0)
            gid_io.WriteNodalResults(
                IS_FLUID,
                (self.model_part).Nodes,
                time,
                0)
            gid_io.WriteNodalResults(
                IS_WATER,
                (self.model_part).Nodes,
                time,
                0)
            gid_io.WriteNodalResults(
                NODAL_H,
                (self.model_part).Nodes,
                time,
                0)
            gid_io.WriteNodalResults(
                DISTANCE,
                (self.model_part).Nodes,
                time,
                0)
            gid_io.WriteNodalResults(
                DISPLACEMENT,
                (self.model_part).Nodes,
                time,
                0)
            gid_io.WriteNodalResults(
                IS_VISITED,
                (self.model_part).Nodes,
                time,
                0)
            gid_io.WriteNodalResults(
                AUX_INDEX,
                (self.model_part).Nodes,
                time,
                0)
            gid_io.PrintOnGaussPoints(IS_WATER_ELEMENT, self.model_part, time)

            gid_io.Flush()
            gid_io.FinalizeResults()

    #
    #
    def CalculateDistanceAndDiviedSet(self, domain_size):

        (self.neigh_finder).ClearNeighbours()
        (self.neigh_finder).Execute()
        distance_tools = ElemBasedDistanceUtilities(self.model_part)
        distance_calculator = BodyDistanceCalculationUtils()

        # assign IS_VISITED1 to elem with DISTANCE>=0 and change DSITANCE to posetive for external ones
        # Assign Zero distance to interface nodes
        for node in (self.model_part).Nodes:
            if(node.GetSolutionStepValue(IS_INTERFACE) == 1.0):
                node.SetSolutionStepValue(DISTANCE, 0, 0.0)

        distance_tools.MarkExternalAndMixedNodes()
        distance_tools.ChangeSignToDistance()

        # calculate distances towards the interior of the domain
        if(domain_size == 2):
            distance_calculator.CalculateDistances2D(
                (self.model_part).Elements,
                DISTANCE,
                True)
        else:
            distance_calculator.CalculateDistances3D(
                (self.model_part).Elements,
                DISTANCE,
                True)

        # change sign
        distance_tools.ChangeSignToDistance()

        # mark as visited all of the nodes inside the fluid domain
        distance_tools.MarkInternalAndMixedNodes()
        print(((self.model_part).Elements).Size())
        # calculate distances towards the outside
        if(domain_size == 2):
            distance_calculator.CalculateDistances2D(
                (self.model_part).Elements,
                DISTANCE,
                True)
        else:
            distance_calculator.CalculateDistances3D(
                (self.model_part).Elements,
                DISTANCE,
                True)

        # Decide IS_WATER flag due to DISTANCE
# for node in (self.model_part).Nodes:
# if(node.GetSolutionStepValue(DISTANCE)<= 0.0):
# node.SetSolutionStepValue(IS_WATER,0,0.0)
# else:
# node.SetSolutionStepValue(IS_WATER,0,1.0)
# if(node.GetSolutionStepValue(DISTANCE)== 0.0):
# print"This node has distance zero, is_interface is assigned"
# node.SetSolutionStepValue(IS_INTERFACE,0,1.0)
# node.SetSolutionStepValue(IS_VISITED,0,1.0)
        # save as distance of the old time step
        distance_tools.SaveScalarVariableToOldStep(DISTANCE)
        print("finished RecalculateDistanceFunction")
     #   (self.SetDivided).SetDividedElem_2D()

        print(">>>>>ELEMENTS ARE DIVIDED<<<<<<<<<<<<")
     #
     #

    def DistToH(self):
        possible_h = self.CalculateRadius()
        print(possible_h)

        min_H = possible_h * 3.14 / 60
        # min_H = .0007#0.001
        sec_min_H = 10 * min_H  # .004
        max_H = .02
        ref_dist = 10 * min_H
        sec_ref_dist = 50 * min_H
        third_ref_dist = 200 * min_H

        slope = (sec_min_H - min_H) / (sec_ref_dist - ref_dist)

        second_slope = (max_H - sec_min_H) / (third_ref_dist - sec_ref_dist)

        # search for min an max of H
# for node in (self.model_part).Nodes:
# node_H = node.GetSolutionStepValue(NODAL_H,0)
# if(node_H<self.min_H):
# self.min_H = node_H
# else:
# if(node_H > self.max_H):
# self.max_H = node_H

        # H = H + dist * dist
        # print ">>>>>DISt TO H ASSIGNMENT<<<<<<<<<<<<"
        for node in (self.model_part).Nodes:
            current_dist = node.GetSolutionStepValue(DISTANCE, 0)
            if(abs(current_dist) <= ref_dist):
                node_H = min_H  # + slope*abs(current_dist)
                node.SetSolutionStepValue(NODAL_H, 0, node_H)

        #    if(ref_dist<abs(current_dist) and abs(current_dist)<= sec_ref_dist):
         #       node_H = min_H + slope*(abs(current_dist) - ref_dist)
          #      node.SetSolutionStepValue(NODAL_H,0,node_H)

        #    if(sec_ref_dist<abs(current_dist) and abs(current_dist)<=third_ref_dist):
        #        node_H = sec_min_H + second_slope*(abs(current_dist)- sec_ref_dist)
         #       node.SetSolutionStepValue(NODAL_H,0,node_H)

        #    if(abs(current_dist)>third_ref_dist):
         #       node_H = max_H
          #      node.SetSolutionStepValue(NODAL_H,0,node_H)

#
    def CalculateRadius(self):
        max_radi = 0.0
        for node in (self.model_part).Nodes:
            if node.GetSolutionStepValue(IS_INTERFACE) == 1.0:
                X_ref = node.X
                Y_ref = node.Y
                Z_ref = node.Z

        for node in (self.model_part).Nodes:
            if node.GetSolutionStepValue(IS_INTERFACE) == 1.0:
                radi = pow(
                    node.X - X_ref,
                    2) + pow(
                        node.Y - Y_ref,
                        2) + pow(
                            node.Z - Z_ref,
                            2)
                if(radi > max_radi):
                    max_radi = radi

        max_radi = pow(max_radi, 0.5)

        if (max_radi == 0.0):
            max_radi = 0.076
        return max_radi



     #
