from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# importing the Kratos Library
import KratosMultiphysics
import KratosMultiphysics.ALEApplication as KratosALE
KratosMultiphysics.CheckForPreviousImport()


def CreateSolver(model_part, custom_settings):
    return MeshSolverBase(model_part, custom_settings)


class MeshSolverBase(object):
    """The base class for mesh motion solvers.
    
    This class defines the user interface to mesh motion solvers.

    Derived classes must override the function _create_mesh_motion_solver()
    to customize the mesh motion algorithm. The mesh motion solver and linear
    solver should always be retrieved using the getter functions. Only the 
    member variables listed below should be accessed directly.

    Public member variables:
    settings -- Kratos parameters for general mesh motion settings.
    model_part -- the mesh motion model part.
    """
    def __init__(self, model_part, custom_settings):
        default_settings = KratosMultiphysics.Parameters("""
        {
            "solver_type" : "mesh_solver_structural_similarity",
            "model_import_settings" : {
                "input_type"     : "mdpa",
                "input_filename" : "unknown_name"
            },
            "ale_linear_solver_settings" : {
                "solver_type" : "AMGCL",
                "smoother_type":"ilu0",
                "krylov_type": "gmres",
                "coarsening_type": "aggregation",
                "max_iteration": 200,
                "provide_coordinates": false,
                "gmres_krylov_space_dimension": 100,
                "verbosity" : 0,
                "tolerance": 1e-7,
                "scaling": false,
                "block_size": 1,
                "use_block_matrices_if_possible" : true,
                "coarse_enough" : 5000
            },
            "time_order" : 2,
            "reform_dofs_each_step" : false,
            "compute_reactions"     : false
        }""")
        self.settings = custom_settings
        self.settings.ValidateAndAssignDefaults(default_settings)
        self.model_part = model_part
        print("::[MeshSolverBase]:: Construction finished")

    #### Public user interface functions ####

    def AddVariables(self):
        self.model_part.AddNodalSolutionStepVariable(KratosMultiphysics.MESH_DISPLACEMENT)
        self.model_part.AddNodalSolutionStepVariable(KratosMultiphysics.MESH_VELOCITY)
        self.model_part.AddNodalSolutionStepVariable(KratosMultiphysics.MESH_REACTION)
        print("::[MeshSolverBase]:: Variables ADDED.")

    def AddDofs(self):
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.MESH_DISPLACEMENT_X, self.model_part)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.MESH_DISPLACEMENT_Y, self.model_part)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.MESH_DISPLACEMENT_Z, self.model_part)
        print("::[MeshSolverBase]:: DOFs ADDED.")

    def Initialize(self):
        self.get_mesh_motion_solver().Initialize()
        print("::[MeshSolverBase]:: Finished initialization.")

    def InitializeSolutionStep(self):
        self.get_mesh_motion_solver().InitializeSolutionStep()

    def FinalizeSolutionStep(self):
        self.get_mesh_motion_solver().FinalizeSolutionStep()

    def SetEchoLevel(self, level):
        self.get_mesh_motion_solver().SetEchoLevel(level)

    def GetEchoLevel(self):
        self.get_mesh_motion_solver().GetEchoLevel()

    def Solve(self):
        self.get_mesh_motion_solver().Solve()

    def Clear(self):
        self.get_mesh_motion_solver().Clear()

    def Check(self):
        self.get_mesh_motion_solver().Check()

    def MoveMesh(self):
        self.get_mesh_motion_solver().MoveMesh()

    def ImportModelPart(self):
        raise Exception("::[MeshSolverBase]:: ImportModelPart() is not implemented yet.")

    def GetComputingModelPart(self):
        return self.model_part

    #### Specific internal functions ####

    def get_linear_solver(self):
        if not hasattr(self, '_linear_solver'):
            self._linear_solver = self._create_linear_solver()
        return self._linear_solver

    def get_mesh_motion_solver(self):
        if not hasattr(self, '_mesh_motion_solver'):
            self._mesh_motion_solver = self._create_mesh_motion_solver()
        return self._mesh_motion_solver

    #### Private functions ####
    
    def _create_linear_solver(self):
        import linear_solver_factory
        linear_solver = linear_solver_factory.ConstructSolver(self.settings["ale_linear_solver_settings"])
        return linear_solver

    def _create_mesh_motion_solver(self):
        """Create the mesh motion solver.
        
        The mesh motion solver must provide the functions defined in SolutionStrategy.
        """
        raise Exception("Mesh motion solver must be created by the derived class.")
