from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# importing the Kratos Library
from KratosMultiphysics import *
from KratosMultiphysics.ALEApplication import *
CheckForPreviousImport()

# import mesh solver base class
import mesh_solver_base

def CreateSolver(model_part, custom_settings):
    return MeshSolverBallvertex(model_part,custom_settings)


class MeshSolverBallvertex(mesh_solver_base.MeshSolverBase):

    def __init__(self, model_part, custom_settings):

        # default settings for ballvertex mesh solver
        default_settings = Parameters("""
        {
            "mesh_reform_dofs_each_step": false
        }""")

        custom_settings.ValidateAndAssignDefaults(default_settings)

        # assign parameters
        self.model_part = model_part
        self.domain_size = model_part.ProcessInfo[DOMAIN_SIZE]
        self.mesh_reform_dofs_each_step = custom_settings["mesh_reform_dofs_each_step"].GetBool()

        # neighbour search
        number_of_avg_elems = 10
        number_of_avg_nodes = 10
        self.neighbour_search = FindNodalNeighboursProcess(model_part, number_of_avg_elems, number_of_avg_nodes)

        # assignation of parameters to be used
        self.time_order = 2

        # definition of the solvers
        # pILUPrecond = ILU0Preconditioner()
        # self.linear_solver =  BICGSTABSolver(1e-5, 300,pILUPrecond)
        # pDiagPrecond = DiagonalPreconditioner()
        # self.linear_solver = CGSolver(1e-3, 300, pDiagPrecond)
        self.linear_solver = ScalingSolver(DeflatedCGSolver(1e-6, 3000, True, 1000), True)

    def AddVariables(self):

        self.model_part.AddNodalSolutionStepVariable(MESH_DISPLACEMENT)
        self.model_part.AddNodalSolutionStepVariable(MESH_VELOCITY)

        print("Mesh solver variables added correctly.")

    def AddDofs(self):

        for node in self.model_part.Nodes:

            node.AddDof(DISPLACEMENT_X)
            node.AddDof(DISPLACEMENT_Y)
            node.AddDof(DISPLACEMENT_Z)

        print("Mesh solver DOFs added correctly.")

    def Initialize(self):
        (self.neighbour_search).Execute()

        if(self.domain_size == 2):
            self.solver = BallVertexMeshMoving2D()
        else:
            self.solver = BallVertexMeshMoving3D()

        self.move_mesh_utilities = MoveMeshUtilities()

        if(self.mesh_reform_dofs_each_step == False):
            self.solver.ConstructSystem(self.model_part)

    def Solve(self):
        if(self.mesh_reform_dofs_each_step):
            (self.neighbour_search).Execute()

            self.solver.ConstructSystem(self.model_part)

            self.solver.BuildAndSolveSystem(self.model_part, self.linear_solver);

            self.solver.ClearSystem()
        else:
            self.solver.BuildAndSolveSystem(self.model_part, self.linear_solver);

        # move mesh and calculate mesh velocity
        self.move_mesh_utilities.BDF_MoveMesh(self.time_order, self.model_part)

    def MoveNodes(self):
        self.move_mesh_utilities.BDF_MoveMesh(self.time_order, self.model_part)
