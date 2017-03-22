from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# importing the Kratos Library
from KratosMultiphysics import *
from KratosMultiphysics.ALEApplication import *
from KratosMultiphysics.ExternalSolversApplication import *
CheckForPreviousImport()

# import mesh solver base class
import mesh_solver_base

def CreateSolver(model_part, custom_settings):
    return MeshSolverLaplacian(model_part, custom_settings)


class MeshSolverLaplacian(mesh_solver_base.MeshSolverBase):

    def __init__(self, model_part, custom_settings):

        # default settings for laplacian mesh solver
        default_settings = Parameters("""
        {
            "time_order"                 : 2,
            "mesh_reform_dofs_each_step" : false
        }""")

        custom_settings.ValidateAndAssignDefaults(default_settings)

        # assign parameters
        self.model_part = model_part
        self.time_order = custom_settings["time_order"].GetInt()
        self.mesh_reform_dofs_each_step = custom_settings["mesh_reform_dofs_each_step"].GetBool()

        # neighbour search
        number_of_avg_elems = 10
        number_of_avg_nodes = 10
        self.neighbour_search = FindNodalNeighboursProcess(model_part, number_of_avg_elems, number_of_avg_nodes)

        # definition of the solvers
        tol = 1e-8
        max_it = 1000
        verbosity = 1
        m = 100
        self.linear_solver = AMGCLSolver(AMGCLSmoother.DAMPED_JACOBI, AMGCLIterativeSolverType.BICGSTAB, tol, max_it, verbosity, m)
        #pILUPrecond = ILU0Preconditioner()
        #self.linear_solver = DeflatedCGSolver(1e-6, 3000, True, 1000)
        #self.linear_solver =  BICGSTABSolver(1e-3, 300,pILUPrecond)
        #self.linear_solver = ScalingSolver(DeflatedCGSolver(1e-6, 3000, True, 1000), True)

    def Initialize(self):
        (self.neighbour_search).Execute()

        self.solver = LaplacianMeshMovingStrategy(self.model_part, self.linear_solver, self.time_order, self.mesh_reform_dofs_each_step)

        (self.solver).SetEchoLevel(0)
        
