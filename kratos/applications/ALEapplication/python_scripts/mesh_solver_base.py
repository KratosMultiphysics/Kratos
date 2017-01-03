from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# importing the Kratos Library
import KratosMultiphysics
import KratosMultiphysics.ALEApplication as KratosALE
KratosMultiphysics.CheckForPreviousImport()


def CreateSolver(model_part, custom_settings):
    return MeshSolverBase(model_part, custom_settings)


class MeshSolverBase:

    def __init__(self, model_part, custom_settings):
        raise Exception("Base mesh solver Python constructor taken. Please implement the __init__ method in your mesh solver.")

    def AddVariables(self):
        self.model_part.AddNodalSolutionStepVariable(KratosALE.MESH_DISPLACEMENT)
        self.model_part.AddNodalSolutionStepVariable(KratosMultiphysics.MESH_VELOCITY)
        self.model_part.AddNodalSolutionStepVariable(KratosALE.MESH_REACTION)
        self.model_part.AddNodalSolutionStepVariable(KratosALE.MESH_RHS)

        print("Mesh base solver variables added correctly.")

    def AddDofs(self):
        for node in self.model_part.Nodes:
            node.AddDof(KratosALE.MESH_DISPLACEMENT_X, KratosALE.MESH_REACTION_X)
            node.AddDof(KratosALE.MESH_DISPLACEMENT_Y, KratosALE.MESH_REACTION_Y)
            node.AddDof(KratosALE.MESH_DISPLACEMENT_Z, KratosALE.MESH_REACTION_Z)

        print("Mesh base solver DOFs added correctly.")

    def Initialize(self):
        raise Exception("Base mesh solver method taken. Please implement the Initialize() method in your mesh solver.")

    def Solve(self):
        if(self.mesh_reform_dofs_each_step):
            (self.neighbour_search).Execute()

        (self.solver).Solve()

    def MoveNodes(self):
        (self.solver).MoveNodes()
