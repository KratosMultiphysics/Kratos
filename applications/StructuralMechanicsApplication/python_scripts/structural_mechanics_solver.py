from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics

# Import applications
import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication

# Importing the base class
from KratosMultiphysics.python_solver import PythonSolver

# Other imports
from KratosMultiphysics.StructuralMechanicsApplication import check_and_prepare_model_process_structural
from KratosMultiphysics.StructuralMechanicsApplication import convergence_criteria_factory
from KratosMultiphysics import python_linear_solver_factory as linear_solver_factory
from KratosMultiphysics import auxiliary_solver_utilities
import KratosMultiphysics.kratos_utilities as kratos_utils

class MechanicalSolver(PythonSolver):
    """The base class for structural mechanics solvers.

    This class provides functions for importing and exporting models,
    adding nodal variables and dofs and solving each solution step.

    Derived classes must override the function _create_solution_scheme which
    constructs and returns a solution scheme. Depending on the type of
    solver, derived classes may also need to override the following functions:

    _create_solution_scheme
    _create_convergence_criterion
    _create_linear_solver
    _create_builder_and_solver
    _create_mechanical_solution_strategy

    The mechanical_solution_strategy, builder_and_solver, etc. should alway be retrieved
    using the getter functions get_mechanical_solution_strategy, get_builder_and_solver,
    etc. from this base class.

    Only the member variables listed below should be accessed directly.

    Public member variables:
    model -- the model containing the modelpart used to construct the solver.
    settings -- Kratos parameters containing solver settings.
    """
    def __init__(self, model, custom_settings):
        settings_have_smps_for_comp_mp = custom_settings.Has("problem_domain_sub_model_part_list") or custom_settings.Has("processes_sub_model_part_list")

        self._validate_settings_in_baseclass=True # To be removed eventually
        super(MechanicalSolver, self).__init__(model, custom_settings)

        model_part_name = self.settings["model_part_name"].GetString()

        if model_part_name == "":
            raise Exception('Please specify a model_part name!')

        # for explicitly constructing the computing modelpart as a submodelpart of the mainmodelpart
        self.use_computing_model_part = custom_settings["use_computing_model_part"].GetBool()
        if not self.use_computing_model_part and settings_have_smps_for_comp_mp:
            raise Exception('"problem_domain_sub_model_part_list" and "processes_sub_model_part_list" can only be specified when NOT using a ComputingModelPart! It is recommended Not to use a ComputingModelPart, then the entire Modelpart is used for the computation. At some point always the entire Modelpart will be used!')

        # Only needed during the transition of removing the ComputingModelPart
        if self.settings["problem_domain_sub_model_part_list"].size() == 0:
            self.settings["problem_domain_sub_model_part_list"].Append(model_part_name)
        if self.settings["processes_sub_model_part_list"].size() == 0:
            self.settings["processes_sub_model_part_list"].Append(model_part_name)

        if self.model.HasModelPart(model_part_name):
            self.main_model_part = self.model[model_part_name]
        else:
            self.main_model_part = self.model.CreateModelPart(model_part_name)
            domain_size = self.settings["domain_size"].GetInt()
            if domain_size < 0:
                raise Exception('Please specify a "domain_size" >= 0!')
            self.main_model_part.ProcessInfo.SetValue(KratosMultiphysics.DOMAIN_SIZE, domain_size)

        KratosMultiphysics.Logger.PrintInfo("::[MechanicalSolver]:: ", "Construction finished")

        # Set if the analysis is restarted
        if self.settings["model_import_settings"]["input_type"].GetString() == "rest":
            self.main_model_part.ProcessInfo[KratosMultiphysics.IS_RESTARTED] = True
        else:
            self.main_model_part.ProcessInfo[KratosMultiphysics.IS_RESTARTED] = False

    @classmethod
    def GetDefaultSettings(cls):
        this_defaults = KratosMultiphysics.Parameters("""{
            "solver_type" : "mechanical_solver",
            "model_part_name" : "",
            "domain_size" : -1,
            "echo_level": 0,
            "buffer_size": 2,
            "analysis_type": "non_linear",
            "model_import_settings": {
                "input_type": "mdpa",
                "input_filename": "unknown_name"
            },
            "computing_model_part_name" : "computing_domain",
            "use_computing_model_part" : true,
            "material_import_settings" :{
                "materials_filename": ""
            },
            "time_stepping" : { },
            "rotation_dofs": false,
            "reform_dofs_at_each_step": false,
            "line_search": false,
            "compute_reactions": true,
            "block_builder": true,
            "clear_storage": false,
            "move_mesh_flag": true,
            "multi_point_constraints_used": true,
            "convergence_criterion": "residual_criterion",
            "displacement_relative_tolerance": 1.0e-4,
            "displacement_absolute_tolerance": 1.0e-9,
            "residual_relative_tolerance": 1.0e-4,
            "residual_absolute_tolerance": 1.0e-9,
            "max_iteration": 10,
            "linear_solver_settings": { },
            "problem_domain_sub_model_part_list": [],
            "processes_sub_model_part_list": [],
            "auxiliary_variables_list" : [],
            "auxiliary_dofs_list" : [],
            "auxiliary_reaction_list" : []
        }""")
        this_defaults.AddMissingParameters(super(MechanicalSolver, cls).GetDefaultSettings())
        return this_defaults

    def AddVariables(self):
        # this can safely be called also for restarts, it is internally checked if the variables exist already
        # Add displacements.
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.REACTION)
        # Add specific variables for the problem conditions.
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.POSITIVE_FACE_PRESSURE)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NEGATIVE_FACE_PRESSURE)
        self.main_model_part.AddNodalSolutionStepVariable(StructuralMechanicsApplication.POINT_LOAD)
        self.main_model_part.AddNodalSolutionStepVariable(StructuralMechanicsApplication.LINE_LOAD)
        self.main_model_part.AddNodalSolutionStepVariable(StructuralMechanicsApplication.SURFACE_LOAD)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VOLUME_ACCELERATION)
        if self.settings["rotation_dofs"].GetBool():
            # Add specific variables for the problem (rotation dofs).
            self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.ROTATION)
            self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.REACTION_MOMENT)
            self.main_model_part.AddNodalSolutionStepVariable(StructuralMechanicsApplication.POINT_MOMENT)
        # Add variables that the user defined in the ProjectParameters
        auxiliary_solver_utilities.AddVariables(self.main_model_part, self.settings["auxiliary_variables_list"])
        KratosMultiphysics.Logger.PrintInfo("::[MechanicalSolver]:: ", "Variables ADDED")

    def GetMinimumBufferSize(self):
        return 2

    def AddDofs(self):
        # this can safely be called also for restarts, it is internally checked if the dofs exist already
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_X, KratosMultiphysics.REACTION_X,self.main_model_part)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_Y, KratosMultiphysics.REACTION_Y,self.main_model_part)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_Z, KratosMultiphysics.REACTION_Z,self.main_model_part)
        if self.settings["rotation_dofs"].GetBool():
            KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.ROTATION_X, KratosMultiphysics.REACTION_MOMENT_X,self.main_model_part)
            KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.ROTATION_Y, KratosMultiphysics.REACTION_MOMENT_Y,self.main_model_part)
            KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.ROTATION_Z, KratosMultiphysics.REACTION_MOMENT_Z,self.main_model_part)

        # Add dofs that the user defined in the ProjectParameters
        auxiliary_solver_utilities.AddDofs(self.main_model_part, self.settings["auxiliary_dofs_list"], self.settings["auxiliary_reaction_list"])
        KratosMultiphysics.Logger.PrintInfo("::[MechanicalSolver]:: ", "DOF's ADDED")

    def ImportModelPart(self):
        """This function imports the ModelPart
        """
        self._ImportModelPart(self.main_model_part, self.settings["model_import_settings"])

    def PrepareModelPart(self):
        if not self.is_restarted():
            # Check and prepare computing model part and import constitutive laws.
            self._execute_after_reading()
            self._set_and_fill_buffer()

        KratosMultiphysics.Logger.PrintInfo("::[MechanicalSolver]::", "ModelPart prepared for Solver.")

    def Initialize(self):
        """Perform initialization after adding nodal variables and dofs to the main model part. """
        KratosMultiphysics.Logger.PrintInfo("::[MechanicalSolver]:: ", "Initializing ...")
        # The mechanical solution strategy is created here if it does not already exist.
        if self.settings["clear_storage"].GetBool():
            self.Clear()
        mechanical_solution_strategy = self.get_mechanical_solution_strategy()
        mechanical_solution_strategy.SetEchoLevel(self.settings["echo_level"].GetInt())
        if not self.is_restarted():
            mechanical_solution_strategy.Initialize()
        else:
            # SetInitializePerformedFlag is not a member of SolvingStrategy but
            # is used by ResidualBasedNewtonRaphsonStrategy.
            try:
                mechanical_solution_strategy.SetInitializePerformedFlag(True)
            except AttributeError:
                pass
        KratosMultiphysics.Logger.PrintInfo("::[MechanicalSolver]:: ", "Finished initialization.")

    def InitializeSolutionStep(self):
        if self.settings["clear_storage"].GetBool():
            self.Clear()
            self.Initialize() #required after clearing
        self.get_mechanical_solution_strategy().InitializeSolutionStep()

    def Predict(self):
        self.get_mechanical_solution_strategy().Predict()

    def SolveSolutionStep(self):
        is_converged = self.get_mechanical_solution_strategy().SolveSolutionStep()
        if not is_converged:
            msg  = "Solver did not converge for step " + str(self.main_model_part.ProcessInfo[KratosMultiphysics.STEP]) + "\n"
            msg += "corresponding to time " + str(self.main_model_part.ProcessInfo[KratosMultiphysics.TIME]) + "\n"
            KratosMultiphysics.Logger.PrintWarning("::[MechanicalSolver]:: ",msg)
        return is_converged

    def FinalizeSolutionStep(self):
        self.get_mechanical_solution_strategy().FinalizeSolutionStep()

    def AdvanceInTime(self, current_time):
        dt = self.ComputeDeltaTime()
        new_time = current_time + dt
        self.main_model_part.ProcessInfo[KratosMultiphysics.STEP] += 1
        self.main_model_part.CloneTimeStep(new_time)

        return new_time

    def ComputeDeltaTime(self):
        return self.settings["time_stepping"]["time_step"].GetDouble()

    def GetComputingModelPart(self):
        if self.use_computing_model_part:
            if not self.main_model_part.HasSubModelPart(self.settings["computing_model_part_name"].GetString()):
                raise Exception("The ComputingModelPart was not created yet!")
            return self.main_model_part.GetSubModelPart(self.settings["computing_model_part_name"].GetString())
        else:
            return self.main_model_part

    def ExportModelPart(self):
        name_out_file = self.settings["model_import_settings"]["input_filename"].GetString()+".out"
        file = open(name_out_file + ".mdpa","w")
        file.close()
        KratosMultiphysics.ModelPartIO(name_out_file, KratosMultiphysics.IO.WRITE).WriteModelPart(self.main_model_part)

    def SetEchoLevel(self, level):
        self.get_mechanical_solution_strategy().SetEchoLevel(level)

    def Clear(self):
        self.get_mechanical_solution_strategy().Clear()

    def Check(self):
        self.get_mechanical_solution_strategy().Check()

    #### Specific internal functions ####

    def get_solution_scheme(self):
        if not hasattr(self, '_solution_scheme'):
            self._solution_scheme = self._create_solution_scheme()
        return self._solution_scheme

    def get_convergence_criterion(self):
        if not hasattr(self, '_convergence_criterion'):
            self._convergence_criterion = self._create_convergence_criterion()
        return self._convergence_criterion

    def get_linear_solver(self):
        if not hasattr(self, '_linear_solver'):
            self._linear_solver = self._create_linear_solver()
        return self._linear_solver

    def get_builder_and_solver(self):
        if (self.settings["multi_point_constraints_used"].GetBool() is False and
            self.GetComputingModelPart().NumberOfMasterSlaveConstraints() > 0):
            self.settings["multi_point_constraints_used"].SetBool(True)
            self._builder_and_solver = self._create_builder_and_solver()
        if not hasattr(self, '_builder_and_solver'):
            self._builder_and_solver = self._create_builder_and_solver()
        return self._builder_and_solver

    def get_mechanical_solution_strategy(self):
        if (self.settings["multi_point_constraints_used"].GetBool() is False and
            self.GetComputingModelPart().NumberOfMasterSlaveConstraints() > 0):
            self._mechanical_solution_strategy = self._create_mechanical_solution_strategy()
        if not hasattr(self, '_mechanical_solution_strategy'):
            self._mechanical_solution_strategy = self._create_mechanical_solution_strategy()
        return self._mechanical_solution_strategy

    def import_constitutive_laws(self):
        materials_filename = self.settings["material_import_settings"]["materials_filename"].GetString()
        if (materials_filename != ""):
            # Add constitutive laws and material properties from json file to model parts.
            material_settings = KratosMultiphysics.Parameters("""{"Parameters": {"materials_filename": ""}} """)
            material_settings["Parameters"]["materials_filename"].SetString(materials_filename)
            KratosMultiphysics.ReadMaterialsUtility(material_settings, self.model)
            materials_imported = True
        else:
            materials_imported = False
        return materials_imported

    def is_restarted(self):
        # this function avoids the long call to ProcessInfo and is also safer
        # in case the detection of a restart is changed later
        return self.main_model_part.ProcessInfo[KratosMultiphysics.IS_RESTARTED]

    #### Private functions ####

    def _execute_after_reading(self):
        """Prepare computing model part and import constitutive laws. """
        if self.use_computing_model_part:
            # construct the computing-modelpart
            # Auxiliary parameters object for the CheckAndPepareModelProcess
            params = KratosMultiphysics.Parameters("{}")
            params.AddValue("model_part_name",self.settings["model_part_name"])
            params.AddValue("computing_model_part_name",self.settings["computing_model_part_name"])
            params.AddValue("problem_domain_sub_model_part_list",self.settings["problem_domain_sub_model_part_list"])
            params.AddValue("processes_sub_model_part_list",self.settings["processes_sub_model_part_list"])
            # Assign mesh entities from domain and process sub model parts to the computing model part.
            check_and_prepare_model_process_structural.CheckAndPrepareModelProcess(self.model, params).Execute()

        # Import constitutive laws.
        materials_imported = self.import_constitutive_laws()
        if materials_imported:
            KratosMultiphysics.Logger.PrintInfo("::[MechanicalSolver]:: ", "Constitutive law was successfully imported.")
        else:
            KratosMultiphysics.Logger.PrintInfo("::[MechanicalSolver]:: ", "Constitutive law was not imported.")

    def _set_and_fill_buffer(self):
        """Prepare nodal solution step data containers and time step information. """
        # Set the buffer size for the nodal solution steps data. Existing nodal
        # solution step data may be lost.
        required_buffer_size = self.settings["buffer_size"].GetInt()
        if required_buffer_size < self.GetMinimumBufferSize():
            required_buffer_size = self.GetMinimumBufferSize()
        current_buffer_size = self.main_model_part.GetBufferSize()
        buffer_size = max(current_buffer_size, required_buffer_size)
        self.main_model_part.SetBufferSize(buffer_size)
        # Cycle the buffer. This sets all historical nodal solution step data to
        # the current value and initializes the time stepping in the process info.
        delta_time = self.main_model_part.ProcessInfo[KratosMultiphysics.DELTA_TIME]
        time = self.main_model_part.ProcessInfo[KratosMultiphysics.TIME]
        step =-buffer_size
        time = time - delta_time * buffer_size
        self.main_model_part.ProcessInfo.SetValue(KratosMultiphysics.TIME, time)
        for i in range(0, buffer_size):
            step = step + 1
            time = time + delta_time
            self.main_model_part.ProcessInfo.SetValue(KratosMultiphysics.STEP, step)
            self.main_model_part.CloneTimeStep(time)

    def _add_dynamic_variables(self):
        # For being consistent for Serial and Trilinos
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VELOCITY)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.ACCELERATION)
        if self.settings["rotation_dofs"].GetBool():
            self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.ANGULAR_VELOCITY)
            self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.ANGULAR_ACCELERATION)

    def _add_dynamic_dofs(self):
        # For being consistent for Serial and Trilinos
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.VELOCITY_X,self.main_model_part)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.VELOCITY_Y,self.main_model_part)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.VELOCITY_Z,self.main_model_part)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.ACCELERATION_X,self.main_model_part)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.ACCELERATION_Y,self.main_model_part)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.ACCELERATION_Z,self.main_model_part)
        if(self.settings["rotation_dofs"].GetBool()):
            KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.ANGULAR_VELOCITY_X,self.main_model_part)
            KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.ANGULAR_VELOCITY_Y,self.main_model_part)
            KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.ANGULAR_VELOCITY_Z,self.main_model_part)
            KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.ANGULAR_ACCELERATION_X,self.main_model_part)
            KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.ANGULAR_ACCELERATION_Y,self.main_model_part)
            KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.ANGULAR_ACCELERATION_Z,self.main_model_part)

    def _get_convergence_criterion_settings(self):
        # Create an auxiliary Kratos parameters object to store the convergence settings.
        conv_params = KratosMultiphysics.Parameters("{}")
        conv_params.AddValue("convergence_criterion",self.settings["convergence_criterion"])
        conv_params.AddValue("rotation_dofs",self.settings["rotation_dofs"])
        conv_params.AddValue("echo_level",self.settings["echo_level"])
        conv_params.AddValue("displacement_relative_tolerance",self.settings["displacement_relative_tolerance"])
        conv_params.AddValue("displacement_absolute_tolerance",self.settings["displacement_absolute_tolerance"])
        conv_params.AddValue("residual_relative_tolerance",self.settings["residual_relative_tolerance"])
        conv_params.AddValue("residual_absolute_tolerance",self.settings["residual_absolute_tolerance"])

        return conv_params

    def _create_convergence_criterion(self):
        convergence_criterion = convergence_criteria_factory.convergence_criterion(self._get_convergence_criterion_settings())
        return convergence_criterion.mechanical_convergence_criterion

    def _create_linear_solver(self):
        linear_solver_configuration = self.settings["linear_solver_settings"]
        if linear_solver_configuration.Has("solver_type"): # user specified a linear solver
            return linear_solver_factory.ConstructSolver(linear_solver_configuration)
        else:
            KratosMultiphysics.Logger.PrintInfo('::[MechanicalSolver]:: No linear solver was specified, using fastest available solver')
            return linear_solver_factory.CreateFastestAvailableDirectLinearSolver()

    def _create_builder_and_solver(self):
        linear_solver = self.get_linear_solver()
        if self.settings["block_builder"].GetBool():
            builder_and_solver = KratosMultiphysics.ResidualBasedBlockBuilderAndSolver(linear_solver)
        else:
            if self.settings["multi_point_constraints_used"].GetBool():
                builder_and_solver = KratosMultiphysics.ResidualBasedEliminationBuilderAndSolverWithConstraints(linear_solver)
            else:
                builder_and_solver = KratosMultiphysics.ResidualBasedEliminationBuilderAndSolver(linear_solver)
        return builder_and_solver

    def _create_solution_scheme(self):
        """Create the solution scheme for the structural problem.
        """
        raise Exception("Solution Scheme creation must be implemented in the derived class.")

    def _create_mechanical_solution_strategy(self):
        analysis_type = self.settings["analysis_type"].GetString()
        if analysis_type == "linear":
            mechanical_solution_strategy = self._create_linear_strategy()
        elif analysis_type == "non_linear":
            if(self.settings["line_search"].GetBool() == False):
                mechanical_solution_strategy = self._create_newton_raphson_strategy()
            else:
                mechanical_solution_strategy = self._create_line_search_strategy()
        else:
            err_msg =  "The requested analysis type \"" + analysis_type + "\" is not available!\n"
            err_msg += "Available options are: \"linear\", \"non_linear\""
            raise Exception(err_msg)
        return mechanical_solution_strategy

    def _create_linear_strategy(self):
        computing_model_part = self.GetComputingModelPart()
        mechanical_scheme = self.get_solution_scheme()
        linear_solver = self.get_linear_solver()
        builder_and_solver = self.get_builder_and_solver()
        return KratosMultiphysics.ResidualBasedLinearStrategy(computing_model_part,
                                                              mechanical_scheme,
                                                              linear_solver,
                                                              builder_and_solver,
                                                              self.settings["compute_reactions"].GetBool(),
                                                              self.settings["reform_dofs_at_each_step"].GetBool(),
                                                              False,
                                                              self.settings["move_mesh_flag"].GetBool())

    def _create_newton_raphson_strategy(self):
        computing_model_part = self.GetComputingModelPart()
        mechanical_scheme = self.get_solution_scheme()
        linear_solver = self.get_linear_solver()
        mechanical_convergence_criterion = self.get_convergence_criterion()
        builder_and_solver = self.get_builder_and_solver()
        return KratosMultiphysics.ResidualBasedNewtonRaphsonStrategy(computing_model_part,
                                                                     mechanical_scheme,
                                                                     linear_solver,
                                                                     mechanical_convergence_criterion,
                                                                     builder_and_solver,
                                                                     self.settings["max_iteration"].GetInt(),
                                                                     self.settings["compute_reactions"].GetBool(),
                                                                     self.settings["reform_dofs_at_each_step"].GetBool(),
                                                                     self.settings["move_mesh_flag"].GetBool())

    def _create_line_search_strategy(self):
        computing_model_part = self.GetComputingModelPart()
        mechanical_scheme = self.get_solution_scheme()
        linear_solver = self.get_linear_solver()
        mechanical_convergence_criterion = self.get_convergence_criterion()
        builder_and_solver = self.get_builder_and_solver()
        return KratosMultiphysics.LineSearchStrategy(computing_model_part,
                                                     mechanical_scheme,
                                                     linear_solver,
                                                     mechanical_convergence_criterion,
                                                     builder_and_solver,
                                                     self.settings["max_iteration"].GetInt(),
                                                     self.settings["compute_reactions"].GetBool(),
                                                     self.settings["reform_dofs_at_each_step"].GetBool(),
                                                     self.settings["move_mesh_flag"].GetBool())
