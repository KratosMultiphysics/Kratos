from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# importing the Kratos Library
from KratosMultiphysics import *
from KratosMultiphysics.ALEApplication import *
from KratosMultiphysics.ExternalSolversApplication import *
CheckForPreviousImport()

# import mesh solver base class
import mesh_solver_base

def CreateSolver(model_part, custom_settings):
    return MeshSolverStructuralSimilarity(model_part, custom_settings)


class MeshSolverStructuralSimilarity(mesh_solver_base.MeshSolverBase):

    def __init__(self, model_part, custom_settings):

        # default settings for structural similarity mesh solver
        default_settings = Parameters("""
        {
            "time_order"                 : 2,
            "mesh_reform_dofs_each_step" : false,
            "mesh_compute_reactions"     : false
        }""")

        custom_settings.ValidateAndAssignDefaults(default_settings)

        # set parameters
        self.model_part = model_part
        self.time_order = custom_settings["time_order"].GetInt()
        self.mesh_compute_reactions = custom_settings["mesh_compute_reactions"].GetBool()
        self.mesh_reform_dofs_each_step = custom_settings["mesh_reform_dofs_each_step"].GetBool()

        # neighbour search
        number_of_avg_elems = 10
        number_of_avg_nodes = 10
        self.neighbour_search = FindNodalNeighboursProcess(model_part, number_of_avg_elems, number_of_avg_nodes)

        # definition of the solvers
        tol = 1e-8
        max_it = 1000
        verbosity = 1
        m = 10
        #self.linear_solver = AMGCLSolver(AMGCLSmoother.GAUSS_SEIDEL, AMGCLIterativeSolverType.BICGSTAB, tol, max_it, verbosity, m)
        #pILUPrecond = ILU0Preconditioner()
        #self.linear_solver =  BICGSTABSolver(1e-9, 300, pILUPrecond)
        self.linear_solver = SuperLUSolver()
        print("Construction of MeshSolverStructuralSimilarity finished")

    def AddVariables(self):
        self.model_part.AddNodalSolutionStepVariable(MESH_DISPLACEMENT)
        self.model_part.AddNodalSolutionStepVariable(MESH_VELOCITY)
        self.model_part.AddNodalSolutionStepVariable(MESH_REACTION)
        self.model_part.AddNodalSolutionStepVariable(MESH_RHS)
        print("Mesh solver variables added correctly.")

    def AddDofs(self):
        for node in self.model_part.Nodes:
            node.AddDof(MESH_DISPLACEMENT_X, MESH_REACTION_X)
            node.AddDof(MESH_DISPLACEMENT_Y, MESH_REACTION_Y)
            node.AddDof(MESH_DISPLACEMENT_Z, MESH_REACTION_Z)
        print("Mesh solver DOFs added correctly.")

    def Initialize(self):
        (self.neighbour_search).Execute()

        # Initialize mesh model part
        for node in self.model_part.Nodes:
            zero = Vector(3)
            zero[0] = 0.0
            zero[1] = 0.0
            zero[2] = 0.0
            node.SetSolutionStepValue(MESH_REACTION,0,zero)
            node.SetSolutionStepValue(MESH_DISPLACEMENT,0,zero)
            node.SetSolutionStepValue(MESH_RHS,0,zero)

        self.solver = StructuralMeshMovingStrategy(self.model_part, self.linear_solver, self.time_order, self.mesh_reform_dofs_each_step, self.mesh_compute_reactions)
        (self.solver).SetEchoLevel(0)

    # def Solve(self):
    #     if(self.mesh_reform_dofs_each_step):
    #         (self.neighbour_search).Execute()
    #     (self.solver).Solve()
    #
    # def MoveNodes(self):
    #     (self.solver).MoveNodes()

    def UpdateReferenceMesh(self):
        (self.solver).UpdateReferenceMesh()
