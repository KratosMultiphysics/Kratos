from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics
import KratosMultiphysics.MeshMovingApplication as KMM

# Other imports
from KratosMultiphysics.python_solver import PythonSolver


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
        self._validate_settings_in_baseclass=True # To be removed eventually
        super(MeshSolverBase,self).__init__(model, custom_settings)

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

        # If required, create the time discretization helper
        if self.settings["calculate_mesh_velocity"].GetBool():
            # BDF2 was the default in the MeshSolver-Strategies
            default_settings = KratosMultiphysics.Parameters("""{
                "time_scheme" : "bdf2",
                "alpha_m": 0.0,
                "alpha_f": 0.0
            }""")

            self.settings["mesh_velocity_calculation"].ValidateAndAssignDefaults(default_settings)
            self.__CreateTimeIntegratorHelper()

        KratosMultiphysics.Logger.PrintInfo("::[MeshSolverBase]:: Construction finished")

    @classmethod
    def GetDefaultSettings(cls):
        this_defaults = KratosMultiphysics.Parameters("""{
            "solver_type"           : "mesh_solver_base",
            "buffer_size"           : 1,
            "domain_size"           : -1,
            "model_part_name"       : "",
            "time_stepping"         : { },
            "model_import_settings" : {
                "input_type"     : "mdpa",
                "input_filename" : "unknown_name"
            },
            "linear_solver_settings" : {
                "solver_type" : "amgcl",
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
            "reform_dofs_each_step"     : false,
            "compute_reactions"         : false,
            "poisson_ratio"             : 0.3,
            "calculate_mesh_velocity"   : true,
            "mesh_velocity_calculation" : { },
            "superimpose_mesh_disp_with": [],
            "superimpose_mesh_velocity_with": []
        }""")
        this_defaults.AddMissingParameters(super(MeshSolverBase, cls).GetDefaultSettings())
        return this_defaults

    #### Public user interface functions ####

    def AddVariables(self):
        # Add variables required for the mesh moving calculation
        self.mesh_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.MESH_DISPLACEMENT)
        self.mesh_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.MESH_REACTION)

        # Adding Variables used for computation of Mesh-Velocity
        if self.settings["calculate_mesh_velocity"].GetBool():
            self.mesh_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.MESH_VELOCITY)
            time_scheme = self.settings["mesh_velocity_calculation"]["time_scheme"].GetString()
            if not time_scheme.startswith("bdf"): # bdfx does not need MESH_ACCELERATION
                self.mesh_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.MESH_ACCELERATION)

        KratosMultiphysics.Logger.PrintInfo("::[MeshSolverBase]:: Variables ADDED.")

    def AddDofs(self):
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.MESH_DISPLACEMENT_X, KratosMultiphysics.MESH_REACTION_X, self.mesh_model_part)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.MESH_DISPLACEMENT_Y, KratosMultiphysics.MESH_REACTION_Y, self.mesh_model_part)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.MESH_DISPLACEMENT_Z, KratosMultiphysics.MESH_REACTION_Z, self.mesh_model_part)
        KratosMultiphysics.Logger.PrintInfo("::[MeshSolverBase]:: DOFs ADDED.")

    def AdvanceInTime(self, current_time):
        dt = self.settings["time_stepping"]["time_step"].GetDouble()
        new_time = current_time + dt
        self.mesh_model_part.ProcessInfo[KratosMultiphysics.STEP] += 1
        self.mesh_model_part.CloneTimeStep(new_time)

        return new_time

    def Initialize(self):
        self.get_mesh_motion_solving_strategy().Initialize()
        #self.neighbour_search.Execute()
        KratosMultiphysics.Logger.PrintInfo("::[MeshSolverBase]:: Finished initialization.")

    def InitializeSolutionStep(self):
        self.get_mesh_motion_solving_strategy().InitializeSolutionStep()

    def FinalizeSolutionStep(self):
        self.get_mesh_motion_solving_strategy().FinalizeSolutionStep()

    def Predict(self):
        self.get_mesh_motion_solving_strategy().Predict()

    def SolveSolutionStep(self):
        # Calling Solve bcs this is what is currently implemented in the MeshSolverStrategies
        # explicit bool conversion is only needed bcs "Solve" returns a double
        is_converged = bool(self.get_mesh_motion_solving_strategy().Solve())
        self.MoveMesh()

        # Superimpose the user-defined mesh displacement
        for variable in KratosMultiphysics.kratos_utilities.GenerateVariableListFromInput(self.settings["superimpose_mesh_disp_with"]):
            KMM.SuperImposeMeshDisplacement(variable)

        # Superimpose the user-defined mesh velocity
        for variable in KratosMultiphysics.kratos_utilities.GenerateVariableListFromInput(self.settings["superimpose_mesh_velocity_with"]):
            KMM.SuperImposeMeshVelocity(variable)

        return is_converged

    def SetEchoLevel(self, level):
        self.get_mesh_motion_solving_strategy().SetEchoLevel(level)

    def GetEchoLevel(self):
        self.get_mesh_motion_solving_strategy().GetEchoLevel()

    def Clear(self):
        self.get_mesh_motion_solving_strategy().Clear()

    def GetMinimumBufferSize(self):
        buffer_size = max(self.settings["buffer_size"].GetInt(), self.mesh_model_part.GetBufferSize())
        if self.settings["calculate_mesh_velocity"].GetBool():
            buffer_size = max(buffer_size, KratosMultiphysics.TimeDiscretization.GetMinimumBufferSize(self.time_int_helper))
        return buffer_size

    def MoveMesh(self):
        # move local and ghost nodes
        self.mesh_model_part.GetCommunicator().SynchronizeVariable(KratosMultiphysics.MESH_DISPLACEMENT)
        KMM.MoveMesh(self.mesh_model_part.Nodes)

        # If required, calculate the MESH_VELOCITY.
        if self.settings["calculate_mesh_velocity"].GetBool():
            KMM.CalculateMeshVelocities(self.mesh_model_part, self.time_int_helper)

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
        from KratosMultiphysics.python_linear_solver_factory import ConstructSolver
        return ConstructSolver(self.settings["linear_solver_settings"])

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

    def __CreateTimeIntegratorHelper(self):
        '''Initializing the helper-class for the time-integration
        '''
        mesh_vel_calc_setting = self.settings["mesh_velocity_calculation"]
        time_scheme = mesh_vel_calc_setting["time_scheme"].GetString()

        if time_scheme == "bdf1":
            self.time_int_helper = KratosMultiphysics.TimeDiscretization.BDF1()
        elif time_scheme == "bdf2":
            self.time_int_helper = KratosMultiphysics.TimeDiscretization.BDF2()
        elif time_scheme == "newmark":
            self.time_int_helper = KratosMultiphysics.TimeDiscretization.Newmark()
        elif time_scheme == "bossak":
            if mesh_vel_calc_setting.Has("alpha_m"):
                alpha_m = mesh_vel_calc_setting["alpha_m"].GetDouble()
                self.time_int_helper = KratosMultiphysics.TimeDiscretization.Bossak(alpha_m)
            else:
                self.time_int_helper = KratosMultiphysics.TimeDiscretization.Bossak()
        elif time_scheme == "generalized_alpha":
            alpha_m = mesh_vel_calc_setting["alpha_m"].GetDouble()
            alpha_f = mesh_vel_calc_setting["alpha_f"].GetDouble()
            self.time_int_helper = KratosMultiphysics.TimeDiscretization.GeneralizedAlpha(alpha_m, alpha_f)
        else:
            err_msg =  'The requested time scheme "{}" is not available for the calculation of the mesh velocity!\n'.format(time_scheme)
            err_msg += 'Available options are: "bdf1", "bdf2", "newmark", "bossak", "generalized_alpha"'
            raise Exception(err_msg)
