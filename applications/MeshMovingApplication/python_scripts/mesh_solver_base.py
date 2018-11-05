from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics

# Check that applications were imported in the main script
KratosMultiphysics.CheckRegisteredApplications("MeshMovingApplication")

# Import applications
import KratosMultiphysics.MeshMovingApplication as KratosMeshMoving

# Other imports
from python_solver import PythonSolver


def CreateSolver(mesh_model_part, custom_settings):
    return MeshSolverBase(mesh_model_part, custom_settings)


class MeshSolverBase(PythonSolver):
    """The base class for mesh motion solvers.

    This class defines the user interface to mesh motion solvers.

    Derived classes must override the function _create_mesh_motion_solving_strategy()
    to customize the mesh motion algorithm. The mesh motion solving strategy and linear
    solver should always be retrieved using the getter functions. Only the
    member variables listed below should be accessed directly.

    Public member variables:
    settings -- Kratos parameters for general mesh motion settings.
    mesh_model_part -- the mesh motion model part.
    """
    def __init__(self, model, custom_settings):
        super(MeshSolverBase,self).__init__(model, custom_settings)

        default_settings = KratosMultiphysics.Parameters("""
        {
            "solver_type"           : "mesh_solver_structural_similarity",
            "buffer_size"           : 3,
            "echo_level"            : 0,
            "domain_size"           : -1,
            "model_part_name"       : "",
            "time_stepping"         : { },
            "model_import_settings" : {
                "input_type"     : "mdpa",
                "input_filename" : "unknown_name"
            },
            "mesh_motion_linear_solver_settings" : {
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
            "reform_dofs_each_step"     : false,
            "compute_reactions"         : false,
            "calculate_mesh_velocities" : true
        }""")

        self.settings.ValidateAndAssignDefaults(default_settings)

        # Either retrieve the model part from the model or create a new one
        model_part_name = self.settings["model_part_name"].GetString()

        if model_part_name == "":
            raise Exception('Please provide the model part name as the "model_part_name" (string) parameter!')

        if self.model.HasModelPart(model_part_name):
            self.mesh_model_part = self.model.GetModelPart(model_part_name)
        else:
            self.mesh_model_part = model.CreateModelPart(model_part_name)

        domain_size = self.settings["domain_size"].GetInt()
        if domain_size == -1:
            raise Exception('Please provide the domain size as the "domain_size" (int) parameter!')

        self.mesh_model_part.ProcessInfo.SetValue(KratosMultiphysics.DOMAIN_SIZE, domain_size)

        self.print_on_rank_zero("::[MeshSolverBase]:: Construction finished")

    #### Public user interface functions ####

    def AddVariables(self):
        self.mesh_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.MESH_DISPLACEMENT)
        self.mesh_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.MESH_REACTION)
        if (self.settings["calculate_mesh_velocities"].GetBool() == True):
            self.mesh_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.MESH_VELOCITY)
        self.print_on_rank_zero("::[MeshSolverBase]:: Variables ADDED.")

    def AddDofs(self):
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.MESH_DISPLACEMENT_X, KratosMultiphysics.MESH_REACTION_X, self.mesh_model_part)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.MESH_DISPLACEMENT_Y, KratosMultiphysics.MESH_REACTION_Y, self.mesh_model_part)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.MESH_DISPLACEMENT_Z, KratosMultiphysics.MESH_REACTION_Z, self.mesh_model_part)
        self.print_on_rank_zero("::[MeshSolverBase]:: DOFs ADDED.")

    def AdvanceInTime(self, current_time):
        dt = self.settings["time_stepping"]["time_step"].GetDouble()
        new_time = current_time + dt
        self.mesh_model_part.ProcessInfo[KratosMultiphysics.STEP] += 1
        self.mesh_model_part.CloneTimeStep(new_time)

        return new_time

    def Initialize(self):
        self.get_mesh_motion_solving_strategy().Initialize()
        #self.neighbour_search.Execute()
        self.print_on_rank_zero("::[MeshSolverBase]:: Finished initialization.")

    def InitializeSolutionStep(self):
        self.get_mesh_motion_solving_strategy().InitializeSolutionStep()

    def FinalizeSolutionStep(self):
        self.get_mesh_motion_solving_strategy().FinalizeSolutionStep()

    def Predict(self):
        self.get_mesh_motion_solving_strategy().Predict()

    def SolveSolutionStep(self):
        self.get_mesh_motion_solving_strategy().Solve() # Calling Solve bcs this is what is currently implemented in the MeshSolverStrategies

    def SetEchoLevel(self, level):
        self.get_mesh_motion_solving_strategy().SetEchoLevel(level)

    def GetEchoLevel(self):
        self.get_mesh_motion_solving_strategy().GetEchoLevel()

    def Clear(self):
        self.get_mesh_motion_solving_strategy().Clear()

    def Check(self):
        self.get_mesh_motion_solving_strategy().Check()

    def GetMinimumBufferSize(self):
        time_order = self.settings["time_order"].GetInt()
        if time_order == 1:
            buffer_size = 2
        elif time_order == 2:
            buffer_size = 3
        else:
            raise Exception('"time_order" can only be 1 or 2!')
        return max(buffer_size, self.settings["buffer_size"].GetInt())

    def MoveMesh(self):
        self.get_mesh_motion_solving_strategy().MoveMesh()

    def ImportModelPart(self):
        # we can use the default implementation in the base class
        self._ImportModelPart(self.mesh_model_part,self.settings["model_import_settings"])

    def PrepareModelPart(self):
        if not self.mesh_model_part.ProcessInfo[KratosMultiphysics.IS_RESTARTED]:
            self._set_and_fill_buffer()

    def GetComputingModelPart(self):
        return self.mesh_model_part

    #### Specific internal functions ####

    def get_linear_solver(self):
        if not hasattr(self, '_linear_solver'):
            self._linear_solver = self._create_linear_solver()
        return self._linear_solver

    def get_mesh_motion_solving_strategy(self):
        if not hasattr(self, '_mesh_motion_solving_strategy'):
            self._mesh_motion_solving_strategy = self._create_mesh_motion_solving_strategy()
        return self._mesh_motion_solving_strategy

    #### Private functions ####

    def _create_linear_solver(self):
        import linear_solver_factory
        linear_solver = linear_solver_factory.ConstructSolver(self.settings["mesh_motion_linear_solver_settings"])
        return linear_solver

    def _create_mesh_motion_solving_strategy(self):
        """Create the mesh motion solving strategy.

        The mesh motion solving strategy must provide the functions defined in SolutionStrategy.
        """
        raise Exception("Mesh motion solving strategy must be created by the derived class.")

    def _set_and_fill_buffer(self):
        """Prepare nodal solution step data containers and time step information. """
        # Set the buffer size for the nodal solution steps data. Existing nodal
        # solution step data may be lost.
        buffer_size = self.GetMinimumBufferSize()
        self.mesh_model_part.SetBufferSize(buffer_size)
        # Cycle the buffer. This sets all historical nodal solution step data to
        # the current value and initializes the time stepping in the process info.
        delta_time = self.mesh_model_part.ProcessInfo[KratosMultiphysics.DELTA_TIME]
        time = self.mesh_model_part.ProcessInfo[KratosMultiphysics.TIME]
        step =-buffer_size
        time = time - delta_time * buffer_size
        self.mesh_model_part.ProcessInfo.SetValue(KratosMultiphysics.TIME, time)
        for i in range(0, buffer_size):
            step = step + 1
            time = time + delta_time
            self.mesh_model_part.ProcessInfo.SetValue(KratosMultiphysics.STEP, step)
            self.mesh_model_part.CloneTimeStep(time)
        self.mesh_model_part.ProcessInfo[KratosMultiphysics.IS_RESTARTED] = False
