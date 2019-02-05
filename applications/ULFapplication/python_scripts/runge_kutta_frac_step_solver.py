from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
#CURRENT IMPLEMENTATION IS JUST FOR 2D
from KratosMultiphysics import *
from KratosMultiphysics.MeshingApplication import *
from KratosMultiphysics.ULFApplication import *

CheckForPreviousImport()

def AddVariables(model_part):
    model_part.AddNodalSolutionStepVariable(VELOCITY)
    model_part.AddNodalSolutionStepVariable(ACCELERATION)
    model_part.AddNodalSolutionStepVariable(PRESSURE)
    model_part.AddNodalSolutionStepVariable(NODAL_MASS)
    model_part.AddNodalSolutionStepVariable(IS_FLUID)
    model_part.AddNodalSolutionStepVariable(IS_STRUCTURE)
    model_part.AddNodalSolutionStepVariable(IS_FREE_SURFACE)
    model_part.AddNodalSolutionStepVariable(IS_BOUNDARY)
    model_part.AddNodalSolutionStepVariable(DISPLACEMENT)
    model_part.AddNodalSolutionStepVariable(VISCOSITY)
    model_part.AddNodalSolutionStepVariable(DENSITY)
    model_part.AddNodalSolutionStepVariable(BODY_FORCE)
    model_part.AddNodalSolutionStepVariable(NODAL_H)

    print("variables for the Runge Kutta Frac Step GLS solver added correctly")


def AddDofs(model_part):
    for node in model_part.Nodes:
        # adding dofs
        node.AddDof(VELOCITY_X)
        node.AddDof(VELOCITY_Y)
        node.AddDof(VELOCITY_Z)

        node.AddDof(PRESSURE)

    print("dofs for the the Runge Kutta Frac Step GLS solver added correctly")


class RungeKuttaFracStepSolver:
    #

    def __init__(self, model_part, domain_size, box_corner1, box_corner2):

        self.model_part = model_part

        self.box_corner1=box_corner1
        self.box_corner2=box_corner2
        #neighbour search
        number_of_avg_elems = 10
        number_of_avg_nodes = 10
        self.neighbour_search = FindNodalNeighboursProcess(model_part,number_of_avg_elems,number_of_avg_nodes)

        # self.move_mesh_strategy = 2
        pDiagPrecond = DiagonalPreconditioner()
        # definition of the solvers
        #pILUPrecond = ILU0Preconditioner()
        self.linear_solver =  BICGSTABSolver(1e-6, 5000,pDiagPrecond)
        print("LIN SOLVER", self.linear_solver)


        self.max_iter = 10

        # default settings
        self.echo_level = 1
        self.CalculateReactionFlag = False
        self.ReformDofSetAtEachStep = True
        self.CalculateNormDxFlag = True
        self.domain_size = domain_size
        # self.MoveMeshFlag = True

        if (self.domain_size == 2):
            self.neigh_finder = FindNodalNeighboursProcess(self.model_part, 9, 18)
        if (self.domain_size == 3):
            self.neigh_finder = FindNodalNeighboursProcess(self.model_part, 20, 30)
        # calculate normals
        self.normal_tools = NormalCalculationUtils()


        self.node_erase_process = NodeEraseProcess(self.model_part);
        self.mark_free_surface_process = MarkFreeSurfaceProcess(self.model_part);

        self.mark_outer_nodes_process = MarkOuterNodesProcess(self.model_part);

        if(domain_size == 2):
            self.Mesher = TriGenPFEMModeler()
            self.fluid_neigh_finder = FindNodalNeighboursProcess(self.model_part,9,18)
            self.condition_neigh_finder = FindConditionsNeighboursProcess(self.model_part,2, 10)
            self.elem_neighbor_finder = FindElementalNeighboursProcess(self.model_part, 2, 10)            
     
        self.mark_fluid_process = MarkFluidProcess(self.model_part);
        self.UlfUtils = UlfUtils()

    #
    def Initialize(self):
        (self.neighbour_search).Execute()
        # calculate the normals to the overall domain
        self.normal_tools.CalculateOnSimplex(
            self.model_part.Conditions,
            self.domain_size)
        # for SLIP condition we need to save these Conditions in a list
        # by now SLIP conditions are identified by FLAG_VARIABLE=3.0. this is
        # done in the constructir of the strategy

        # creating the solution strategy
        if (self.domain_size == 2):
            self.solver = RungeKuttaFracStepStrategy2D(
                self.model_part, self.linear_solver, self.CalculateReactionFlag,
                self.ReformDofSetAtEachStep, self.CalculateNormDxFlag)
        if (self.domain_size == 3):
            msg = """Fluid solver error: Your domain size is 3D, but this model is implemented in 2D only to date."""
            raise Exception(msg)

        (self.solver).SetEchoLevel(self.echo_level)
        
        (self.neigh_finder).Execute()
        (self.fluid_neigh_finder).Execute();
        (self.mark_free_surface_process).Execute();
        for node in self.model_part.Nodes:
            node.SetSolutionStepValue(IS_FREE_SURFACE,0,0.0)
        
        Hfinder  = FindNodalHProcess(self.model_part);
        Hfinder.Execute();
        self.Remesh(); 

    def Solve(self):

        (self.solver).SolveStep1()
        #re-meshing is performed after the explicit step and prior to the implicit one (pressure)
        self.Remesh();

        (self.neigh_finder).Execute()
        #setting pressure to zero at free-surface to solve the Poisson's equation        
        for node in (self.model_part).Nodes:
            node.Free(PRESSURE)
            if(node.GetSolutionStepValue(IS_FREE_SURFACE)== 1.0):
                node.SetSolutionStepValue(PRESSURE,0,0.0)
                node.Fix(PRESSURE)


        (self.solver).SolveStep2()

        (self.solver).SolveStep3()

    #
    def SetEchoLevel(self, level):
        (self.solver).SetEchoLevel(level)

    def Remesh(self):	

        alpha_shape=1.4;
        h_factor=0.2

        if(self.domain_size == 2):
            for node in (self.model_part).Nodes: 
                node.SetSolutionStepValue(NODAL_H,0,0.002) 


        self.node_erase_process = NodeEraseProcess(self.model_part);
        self.UlfUtils.MarkLonelyNodesForErasing(self.model_part)

        (self.mark_outer_nodes_process).MarkOuterNodes(self.box_corner1, self.box_corner2);
        self.UlfUtils.MarkNodesCloseToWall(self.model_part, self.domain_size, 2.5000)
        self.UlfUtils.MarkExcessivelyCloseNodes(self.model_part.Nodes, 0.2)
        self.node_erase_process.Execute()

        ((self.model_part).Elements).clear();
        ((self.model_part).Conditions).clear();

        if (self.domain_size == 2):
            (self.Mesher).ReGenerateMesh("Fluid2DGLS_expl","Condition2D", self.model_part, self.node_erase_process, True, True, alpha_shape, h_factor)
             
        for node in (self.model_part).Nodes:
            node.Set(TO_ERASE, False)

        (self.fluid_neigh_finder).Execute();
        (self.elem_neighbor_finder).Execute()
        (self.condition_neigh_finder).Execute();

        (self.mark_free_surface_process).Execute(); 

 

