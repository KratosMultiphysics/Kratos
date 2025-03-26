# Importing the Kratos Library
import KratosMultiphysics

# Import applications
import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication

# Importing the base class
from KratosMultiphysics.python_solver import PythonSolver

# Other imports
from KratosMultiphysics.StructuralMechanicsApplication import convergence_criteria_factory
from KratosMultiphysics import python_linear_solver_factory as linear_solver_factory
from KratosMultiphysics import auxiliary_solver_utilities
from KratosMultiphysics import kratos_utilities

# Other imports
from importlib import import_module

class MechanicalSolver(PythonSolver):
    """The base class for structural mechanics solvers.

    This class provides functions for importing and exporting models,
    adding nodal variables and dofs and solving each solution step.

    Derived classes must override the function _CreateScheme which
    constructs and returns a solution scheme. Depending on the type of
    solver, derived classes may also need to override the following functions:

    _CreateScheme
    _CreateConvergenceCriterion
    _CreateLinearSolver
    _CreateBuilderAndSolver
    _CreateSolutionStrategy

    The mechanical_solution_strategy, builder_and_solver, etc. should alway be retrieved
    using the getter functions _GetSolutionStrategy, get_builder_and_solver,
    etc. from this base class.

    Only the member variables listed below should be accessed directly.

    Public member variables:
    model -- the model containing the modelpart used to construct the solver.
    settings -- Kratos parameters containing solver settings.
    """
    def __init__(self, model, custom_settings):
        self.__ConvertDeprecatedSettings(custom_settings)
        self._validate_settings_in_baseclass=True # To be removed eventually
        super().__init__(model, custom_settings)

        model_part_name = self.settings["model_part_name"].GetString()

        if model_part_name == "":
            raise Exception('Please specify a model_part name!')

        if self.model.HasModelPart(model_part_name):
            self.main_model_part = self.model[model_part_name]
        else:
            self.main_model_part = self.model.CreateModelPart(model_part_name)
            domain_size = self.settings["domain_size"].GetInt()
            if domain_size < 0:
                raise Exception('Please specify a "domain_size" >= 0!')
            self.main_model_part.ProcessInfo.SetValue(KratosMultiphysics.DOMAIN_SIZE, domain_size)

        # Set the Orthogonal SubScales (OSS) switch
        # Note that this needs to be done before the scheme initialize in order to allocate the projection variables
        self.main_model_part.ProcessInfo[KratosMultiphysics.OSS_SWITCH] = self.settings["use_orthogonal_subscales"].GetBool()

        # Some variables initialization
        self.mpc_block_builder_initialized = False

        # Set if the analysis is restarted
        if self.settings["model_import_settings"]["input_type"].GetString() == "rest":
            self.main_model_part.ProcessInfo[KratosMultiphysics.IS_RESTARTED] = True
        else:
            self.main_model_part.ProcessInfo[KratosMultiphysics.IS_RESTARTED] = False

        # Printing message
        KratosMultiphysics.Logger.PrintInfo("::[MechanicalSolver]:: ", "Construction finished")

    @classmethod
    def GetDefaultParameters(cls):
        this_defaults = KratosMultiphysics.Parameters("""{
            "solver_type" : "mechanical_solver",
            "model_part_name" : "",
            "computing_sub_model_part_name" : "",
            "domain_size" : -1,
            "echo_level": 0,
            "buffer_size": 2,
            "analysis_type": "non_linear",
            "model_import_settings": {
                "input_type": "mdpa"
            },
            "material_import_settings" :{
                "materials_filename": ""
            },
            "time_stepping" : { },
            "volumetric_strain_dofs": false,
            "strain_dofs": false,
            "rotation_dofs": false,
            "pressure_dofs": false,
            "displacement_control": false,
            "reform_dofs_at_each_step": false,
            "use_old_stiffness_in_first_iteration": false,
            "use_orthogonal_subscales" : false,
            "compute_reactions": true,
            "solving_strategy_settings": {
                "type" : "newton_raphson",
                "advanced_settings" : { }
            },
            "builder_and_solver_settings" : {
                "type" : "block",
                "advanced_settings" : { }
            },
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
            "auxiliary_variables_list" : [],
            "auxiliary_dofs_list" : [],
            "auxiliary_reaction_list" : []
        }""")
        this_defaults.AddMissingParameters(super().GetDefaultParameters())
        return this_defaults

    def ValidateSettings(self):
        """This function validates the settings of the solver
        """
        super().ValidateSettings()

        # Validate some subparameters
        self.settings["builder_and_solver_settings"].ValidateAndAssignDefaults(self.GetDefaultParameters()["builder_and_solver_settings"])

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
        if self.settings["volumetric_strain_dofs"].GetBool():
            # Add specific variables for the problem (volumetric strain dofs).
            self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VOLUMETRIC_STRAIN)
            self.main_model_part.AddNodalSolutionStepVariable(StructuralMechanicsApplication.REACTION_STRAIN)
            #TODO: These are not required in the standard ASGS case
            #TODO: We can get rid of this overhead once we move to the specification-based variables and DOFs addition
            self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT_PROJECTION)
            self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VOLUMETRIC_STRAIN_PROJECTION)
            self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT_PROJECTION_REACTION)
            self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VOLUMETRIC_STRAIN_PROJECTION_REACTION)
        if self.settings["strain_dofs"].GetBool():
            # Add specific variables for the problem (strain Voigt notation components dofs).
            dim = self.settings["domain_size"].GetInt()
            self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.STRAIN_VECTOR_XX)
            self.main_model_part.AddNodalSolutionStepVariable(StructuralMechanicsApplication.REACTION_STRAIN_VECTOR_XX)
            self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.STRAIN_VECTOR_YY)
            self.main_model_part.AddNodalSolutionStepVariable(StructuralMechanicsApplication.REACTION_STRAIN_VECTOR_YY)
            if dim == 3:
                self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.STRAIN_VECTOR_ZZ)
                self.main_model_part.AddNodalSolutionStepVariable(StructuralMechanicsApplication.REACTION_STRAIN_VECTOR_ZZ)
            self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.STRAIN_VECTOR_XY)
            self.main_model_part.AddNodalSolutionStepVariable(StructuralMechanicsApplication.REACTION_STRAIN_VECTOR_XY)
            if dim == 3:
                self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.STRAIN_VECTOR_YZ)
                self.main_model_part.AddNodalSolutionStepVariable(StructuralMechanicsApplication.REACTION_STRAIN_VECTOR_YZ)
                self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.STRAIN_VECTOR_XZ)
                self.main_model_part.AddNodalSolutionStepVariable(StructuralMechanicsApplication.REACTION_STRAIN_VECTOR_XZ)
        if self.settings["displacement_control"].GetBool():
            # Add displacement-control variables
            self.main_model_part.AddNodalSolutionStepVariable(StructuralMechanicsApplication.LOAD_FACTOR)
            self.main_model_part.AddNodalSolutionStepVariable(StructuralMechanicsApplication.PRESCRIBED_DISPLACEMENT)
        # Add variables that the user defined in the ProjectParameters
        auxiliary_solver_utilities.AddVariables(self.main_model_part, self.settings["auxiliary_variables_list"])
        KratosMultiphysics.Logger.PrintInfo("::[MechanicalSolver]:: ", "Variables ADDED")

    def GetMinimumBufferSize(self):
        return 2

    def AddDofs(self):
        # Append formulation-related DOFs and reactions
        dofs_and_reactions_to_add = []
        dofs_and_reactions_to_add.append(["DISPLACEMENT_X", "REACTION_X"])
        dofs_and_reactions_to_add.append(["DISPLACEMENT_Y", "REACTION_Y"])
        dofs_and_reactions_to_add.append(["DISPLACEMENT_Z", "REACTION_Z"])
        if self.settings["rotation_dofs"].GetBool():
            dofs_and_reactions_to_add.append(["ROTATION_X", "REACTION_MOMENT_X"])
            dofs_and_reactions_to_add.append(["ROTATION_Y", "REACTION_MOMENT_Y"])
            dofs_and_reactions_to_add.append(["ROTATION_Z", "REACTION_MOMENT_Z"])
        if self.settings["volumetric_strain_dofs"].GetBool():
            dofs_and_reactions_to_add.append(["VOLUMETRIC_STRAIN", "REACTION_STRAIN"])
            #TODO: These are only required in the nonlinear OSS case so we are adding them for nothing in the linearised OSS and ASGS
            #TODO: We can get rid of this overhead once we move to the specification-based variables and DOFs addition
            dofs_and_reactions_to_add.append(["DISPLACEMENT_PROJECTION_X", "DISPLACEMENT_PROJECTION_REACTION_X"])
            dofs_and_reactions_to_add.append(["DISPLACEMENT_PROJECTION_Y", "DISPLACEMENT_PROJECTION_REACTION_Y"])
            dofs_and_reactions_to_add.append(["DISPLACEMENT_PROJECTION_Z", "DISPLACEMENT_PROJECTION_REACTION_Z"])
            dofs_and_reactions_to_add.append(["VOLUMETRIC_STRAIN_PROJECTION", "VOLUMETRIC_STRAIN_PROJECTION_REACTION"])
        if self.settings["strain_dofs"].GetBool():
            dim = self.settings["domain_size"].GetInt()
            dofs_and_reactions_to_add.append(["STRAIN_VECTOR_XX", "REACTION_STRAIN_VECTOR_XX"])
            dofs_and_reactions_to_add.append(["STRAIN_VECTOR_YY", "REACTION_STRAIN_VECTOR_YY"])
            if dim == 3:
                dofs_and_reactions_to_add.append(["STRAIN_VECTOR_ZZ", "REACTION_STRAIN_VECTOR_ZZ"])
            dofs_and_reactions_to_add.append(["STRAIN_VECTOR_XY", "REACTION_STRAIN_VECTOR_XY"])
            if dim == 3:
                dofs_and_reactions_to_add.append(["STRAIN_VECTOR_YZ", "REACTION_STRAIN_VECTOR_YZ"])
                dofs_and_reactions_to_add.append(["STRAIN_VECTOR_XZ", "REACTION_STRAIN_VECTOR_XZ"])
        if self.settings["displacement_control"].GetBool():
            dofs_and_reactions_to_add.append(["LOAD_FACTOR", "PRESCRIBED_DISPLACEMENT"])

        # Append user-defined DOFs and reactions in the ProjectParameters
        auxiliary_solver_utilities.AddAuxiliaryDofsToDofsWithReactionsList(
            self.settings["auxiliary_dofs_list"],
            self.settings["auxiliary_reaction_list"],
            dofs_and_reactions_to_add)

        KratosMultiphysics.VariableUtils.AddDofsList(dofs_and_reactions_to_add, self.main_model_part)
        KratosMultiphysics.Logger.PrintInfo("::[MechanicalSolver]:: ", "DOF's ADDED")

    def GetDofsList(self):
        """This function creates and returns a list with the DOFs defined in the conditions and elements specifications
        """
        return KratosMultiphysics.SpecificationsUtilities.GetDofsListFromSpecifications(self.main_model_part)

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
        mechanical_solution_strategy = self._GetSolutionStrategy()
        mechanical_solution_strategy.SetEchoLevel(self.settings["echo_level"].GetInt())
        mechanical_solution_strategy.Initialize()

        # Printing that inialization is finished
        KratosMultiphysics.Logger.PrintInfo("::[MechanicalSolver]:: ", "Finished initialization.")

    def InitializeSolutionStep(self):
        if self.settings["clear_storage"].GetBool():
            self.Clear()
            self.Initialize() #required after clearing
        self._GetSolutionStrategy().InitializeSolutionStep()

    def Predict(self):
        self._GetSolutionStrategy().Predict()

    def SolveSolutionStep(self):
        is_converged = self._GetSolutionStrategy().SolveSolutionStep()
        if not is_converged:
            msg  = "Solver did not converge for step " + str(self.main_model_part.ProcessInfo[KratosMultiphysics.STEP]) + "\n"
            msg += "corresponding to time " + str(self.main_model_part.ProcessInfo[KratosMultiphysics.TIME]) + "\n"
            KratosMultiphysics.Logger.PrintWarning("::[MechanicalSolver]:: ",msg)
        return is_converged

    def FinalizeSolutionStep(self):
        self._GetSolutionStrategy().FinalizeSolutionStep()

    def AdvanceInTime(self, current_time):
        dt = self.ComputeDeltaTime()
        new_time = current_time + dt
        self.main_model_part.ProcessInfo[KratosMultiphysics.STEP] += 1
        self.main_model_part.CloneTimeStep(new_time)

        return new_time

    def ComputeDeltaTime(self):
        if self.settings["time_stepping"].Has("time_step"):
            return self.settings["time_stepping"]["time_step"].GetDouble()
        elif self.settings["time_stepping"].Has("time_step_table"):
            current_time = self.main_model_part.ProcessInfo[KratosMultiphysics.TIME]
            time_step_table = self.settings["time_stepping"]["time_step_table"].GetMatrix()
            tb = KratosMultiphysics.PiecewiseLinearTable()
            for interval in range(time_step_table.Size1()):
                tb.AddRow(time_step_table[interval, 0], time_step_table[interval, 1])
            return tb.GetValue(current_time)
        else:
            raise Exception("::[MechanicalSolver]:: Time stepping not defined!")

    def GetComputingModelPart(self):
        computing_sub_model_part_name = self.settings["computing_sub_model_part_name"].GetString()
        if computing_sub_model_part_name == "":
            # if the user didn't specify a SubModelPart, then use the MainModelPart
            return self.main_model_part
        else:
            computing_model_part_name = self.main_model_part.Name + "." + computing_sub_model_part_name
            return self.model[computing_model_part_name]

    def ExportModelPart(self):
        name_out_file = self.settings["model_import_settings"]["input_filename"].GetString()+".out"
        file = open(name_out_file + ".mdpa","w")
        file.close()
        KratosMultiphysics.ModelPartIO(name_out_file, KratosMultiphysics.IO.WRITE).WriteModelPart(self.main_model_part)

    def SetEchoLevel(self, level):
        self._GetSolutionStrategy().SetEchoLevel(level)

    def Clear(self):
        self._GetSolutionStrategy().Clear()

    def Check(self):
        self._GetSolutionStrategy().Check()

    #### Specific internal functions ####

    def _GetScheme(self):
        if not hasattr(self, '_solution_scheme'):
            self._solution_scheme = self._CreateScheme()
        return self._solution_scheme

    def _GetConvergenceCriterion(self):
        if not hasattr(self, '_convergence_criterion'):
            self._convergence_criterion = self._CreateConvergenceCriterion()
        return self._convergence_criterion

    def _GetLinearSolver(self):
        if not hasattr(self, '_linear_solver'):
            self._linear_solver = self._CreateLinearSolver()
        return self._linear_solver

    def _GetBuilderAndSolver(self):
        builder_and_solver_type: str = self.settings["builder_and_solver_settings"]["type"].GetString()
        if not hasattr(self, '_builder_and_solver'):
            self._builder_and_solver = self._CreateBuilderAndSolver()
        elif builder_and_solver_type == "elimination": # Block builder and solver are unified with MPC and without. In the case of the elimination this could be a problem
            if self.GetComputingModelPart().NumberOfMasterSlaveConstraints() > 0 and not self.mpc_block_builder_initialized:
                self.settings["multi_point_constraints_used"].SetBool(True)
                self._builder_and_solver = self._CreateBuilderAndSolver()
                self.mpc_block_builder_initialized = True
        return self._builder_and_solver

    def _GetSolutionStrategy(self):
        builder_and_solver_type: str = self.settings["builder_and_solver_settings"]["type"].GetString()
        if not hasattr(self, '_mechanical_solution_strategy'):
            self._mechanical_solution_strategy = self._CreateSolutionStrategy()
        elif builder_and_solver_type == "elimination": # Block builder and solver are unified with MPC and without. In the case of the elimination this could be a problem
            if self.GetComputingModelPart().NumberOfMasterSlaveConstraints() > 0 and not self.mpc_block_builder_initialized:
                self._mechanical_solution_strategy = self._CreateSolutionStrategy()
        return self._mechanical_solution_strategy

    def import_constitutive_laws(self):
        if self.settings["material_import_settings"].Has("custom_reader"): # We use our own file for reading
            custom_reader = import_module(self.settings["material_import_settings"]["custom_reader"].GetString())
            custom_reader.ReadMaterials(self.model, self.settings["material_import_settings"])
            materials_imported = True
        else: # We follow the normal path
            materials_filename = self.settings["material_import_settings"]["materials_filename"].GetString()
            if materials_filename != "":
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
        """Import constitutive laws."""
        # Import constitutive laws.
        materials_imported = self.import_constitutive_laws()
        if materials_imported:
            KratosMultiphysics.Logger.PrintInfo("::[MechanicalSolver]:: ", "Constitutive law was successfully imported.")
        else:
            KratosMultiphysics.Logger.PrintInfo("::[MechanicalSolver]:: ", "Constitutive law was not imported.")

    def _set_and_fill_buffer(self):
        required_buffer_size = self.settings["buffer_size"].GetInt()
        if required_buffer_size < self.GetMinimumBufferSize():
            required_buffer_size = self.GetMinimumBufferSize()
        delta_time = self.main_model_part.ProcessInfo[KratosMultiphysics.DELTA_TIME]
        auxiliary_solver_utilities.SetAndFillBuffer(self.main_model_part, required_buffer_size, delta_time)

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
        if self.settings["rotation_dofs"].GetBool():
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
        conv_params.AddValue("volumetric_strain_dofs",self.settings["volumetric_strain_dofs"])
        conv_params.AddValue("strain_dofs",self.settings["strain_dofs"])
        conv_params.AddValue("rotation_dofs",self.settings["rotation_dofs"])
        conv_params.AddValue("echo_level",self.settings["echo_level"])
        conv_params.AddValue("displacement_relative_tolerance",self.settings["displacement_relative_tolerance"])
        conv_params.AddValue("displacement_absolute_tolerance",self.settings["displacement_absolute_tolerance"])
        conv_params.AddValue("residual_relative_tolerance",self.settings["residual_relative_tolerance"])
        conv_params.AddValue("residual_absolute_tolerance",self.settings["residual_absolute_tolerance"])

        return conv_params

    def _CreateConvergenceCriterion(self):
        convergence_criterion = convergence_criteria_factory.convergence_criterion(self._get_convergence_criterion_settings())
        return convergence_criterion.mechanical_convergence_criterion

    def _CreateLinearSolver(self):
        linear_solver_configuration = self.settings["linear_solver_settings"]
        if linear_solver_configuration.Has("solver_type"): # user specified a linear solver
            return linear_solver_factory.ConstructSolver(linear_solver_configuration)
        else:
            KratosMultiphysics.Logger.PrintInfo('::[MechanicalSolver]:: No linear solver was specified, using fastest available solver')
            return linear_solver_factory.CreateFastestAvailableDirectLinearSolver()

    def _CreateBuilderAndSolver(self) -> KratosMultiphysics.BuilderAndSolver:
        if not self.settings["builder_and_solver_settings"].Has("type"):
            return self.__DeprecatedBuilderAndSolverFactory()

        linear_solver: KratosMultiphysics.LinearSolver = self._GetLinearSolver()
        builder_and_solver_name: str = self.settings["builder_and_solver_settings"]["type"].GetString()

        if builder_and_solver_name == "block":
            bs_params = self.settings["builder_and_solver_settings"]["advanced_settings"]
            return KratosMultiphysics.ResidualBasedBlockBuilderAndSolver(linear_solver, bs_params)
        elif builder_and_solver_name == "block_lagrange":
            return KratosMultiphysics.ResidualBasedBlockBuilderAndSolverWithLagrangeMultiplier(linear_solver, bs_params)
        elif builder_and_solver_name == "elimination":
            if self.settings["multi_point_constraints_used"].GetBool():
                return KratosMultiphysics.ResidualBasedEliminationBuilderAndSolverWithConstraints(linear_solver)
            else:
                return KratosMultiphysics.ResidualBasedEliminationBuilderAndSolver(linear_solver)
        else:
            message: str = f"Invalid type for builder and solver '{builder_and_solver_name}'. Options are:\n"
            message += "'block'\n"
            message += "'block_lagrange'\n"
            message += "'elimination'"
            raise ValueError(message)

    def _CreateScheme(self):
        """Create the solution scheme for the structural problem.
        """
        raise Exception("Solution Scheme creation must be implemented in the derived class.")

    def _CreateSolutionStrategy(self):
        analysis_type = self.settings["analysis_type"].GetString()
        if analysis_type == "linear":
            mechanical_solution_strategy = self._create_linear_strategy()
        elif analysis_type == "non_linear":
            # Create strategy
            if self.settings["solving_strategy_settings"]["type"].GetString() == "newton_raphson":
                mechanical_solution_strategy = self._create_newton_raphson_strategy()
            elif self.settings["solving_strategy_settings"]["type"].GetString() == "line_search":
                mechanical_solution_strategy = self._create_line_search_strategy()
            elif self.settings["solving_strategy_settings"]["type"].GetString() == "arc_length":
                mechanical_solution_strategy = self._create_arc_length_strategy()

        else:
            err_msg =  "The requested analysis type \"" + analysis_type + "\" is not available!\n"
            err_msg += "Available options are: \"linear\", \"non_linear\""
            raise Exception(err_msg)
        return mechanical_solution_strategy

    def _create_linear_strategy(self):
        computing_model_part = self.GetComputingModelPart()
        mechanical_scheme = self._GetScheme()
        builder_and_solver = self._GetBuilderAndSolver()
        return KratosMultiphysics.ResidualBasedLinearStrategy(computing_model_part,
                                                              mechanical_scheme,
                                                              builder_and_solver,
                                                              self.settings["compute_reactions"].GetBool(),
                                                              self.settings["reform_dofs_at_each_step"].GetBool(),
                                                              False,
                                                              self.settings["move_mesh_flag"].GetBool())

    def _create_newton_raphson_strategy(self):
        computing_model_part = self.GetComputingModelPart()
        mechanical_scheme = self._GetScheme()
        mechanical_convergence_criterion = self._GetConvergenceCriterion()
        builder_and_solver = self._GetBuilderAndSolver()
        strategy = KratosMultiphysics.ResidualBasedNewtonRaphsonStrategy(computing_model_part,
                                                                     mechanical_scheme,
                                                                     mechanical_convergence_criterion,
                                                                     builder_and_solver,
                                                                     self.settings["max_iteration"].GetInt(),
                                                                     self.settings["compute_reactions"].GetBool(),
                                                                     self.settings["reform_dofs_at_each_step"].GetBool(),
                                                                     self.settings["move_mesh_flag"].GetBool())
        strategy.SetUseOldStiffnessInFirstIterationFlag(self.settings["use_old_stiffness_in_first_iteration"].GetBool())
        return strategy

    def _create_line_search_strategy(self):
        computing_model_part = self.GetComputingModelPart()
        mechanical_scheme = self._GetScheme()
        linear_solver = self._GetLinearSolver()
        mechanical_convergence_criterion = self._GetConvergenceCriterion()
        builder_and_solver = self._GetBuilderAndSolver()
        strategy = KratosMultiphysics.LineSearchStrategy(computing_model_part,
                                                     mechanical_scheme,
                                                     linear_solver,
                                                     mechanical_convergence_criterion,
                                                     builder_and_solver,
                                                     self.settings["max_iteration"].GetInt(),
                                                     self.settings["compute_reactions"].GetBool(),
                                                     self.settings["reform_dofs_at_each_step"].GetBool(),
                                                     self.settings["move_mesh_flag"].GetBool())
        strategy.SetUseOldStiffnessInFirstIterationFlag(self.settings["use_old_stiffness_in_first_iteration"].GetBool())
        return strategy

    def _create_arc_length_strategy(self):
        settings = self.settings["solving_strategy_settings"]["advanced_settings"]
        settings.AddValue("max_iteration", self.settings["max_iteration"])
        settings.AddValue("compute_reactions", self.settings["compute_reactions"])
        settings.AddValue("reform_dofs_at_each_step", self.settings["reform_dofs_at_each_step"])
        settings.AddValue("move_mesh_flag", self.settings["move_mesh_flag"])
        solving_strategy = KratosMultiphysics.ArcLengthStrategy(self.GetComputingModelPart(),
                                                                self._GetScheme(),
                                                                self._GetConvergenceCriterion(),
                                                                self._GetBuilderAndSolver(),
                                                                settings)
        return solving_strategy

    def __ConvertDeprecatedSettings(self, settings: KratosMultiphysics.Parameters) -> None:
        old_unused_settings = [
            "use_computing_model_part",
            "computing_model_part_name",
            "problem_domain_sub_model_part_list",
            "processes_sub_model_part_list"
        ]

        for old_setting in old_unused_settings:
            if settings.Has(old_setting):
                KratosMultiphysics.Logger.PrintWarning("::[MechanicalSolver]:: ", 'Settings contain no longer used setting, please remove it: "{}"'.format(old_setting))
                settings.RemoveValue(old_setting)

        settings_have_use_block_builder = settings.Has("block_builder")
        if settings_have_use_block_builder:
            kratos_utilities.IssueDeprecationWarning('MechanicalSolver', 'Using "block_builder", please move it to "builder_and_solver_settings" as "use_block_builder"')
            if not settings.Has("builder_and_solver_settings"):
                settings.AddEmptyValue("builder_and_solver_settings")

            settings["builder_and_solver_settings"].AddValue("use_block_builder", settings["block_builder"])
            settings.RemoveValue("block_builder")

        settings_have_line_search = settings.Has("line_search")
        if settings_have_line_search:
            kratos_utilities.IssueDeprecationWarning('MechanicalSolver', 'Using "line_search", please move it to "solving_strategy_settings" as "type"')
            if settings["line_search"].GetBool():
                if not settings.Has("solving_strategy_settings"):
                    settings.AddEmptyValue("solving_strategy_settings")

                settings["solving_strategy_settings"].AddEmptyValue("type")
                settings["solving_strategy_settings"]["type"].SetString("line_search")
            settings.RemoveValue("line_search")

        if settings.Has("builder_and_solver_settings") and not settings["builder_and_solver_settings"].Has("type"):
            kratos_utilities.IssueDeprecationWarning(
                "MechanicalSolver",
                "Using deprecated builder and solver settings. Provide 'type' and 'advanced_settings' in the new system.")
            bs_settings = settings["builder_and_solver_settings"]
            updated_bs_settings = KratosMultiphysics.Parameters("""{}""")
            if bs_settings.Has("use_block_builder"):
                if bs_settings["use_block_builder"].GetBool():
                    if bs_settings.Has("use_lagrange_BS") and bs_settings["use_lagrange_BS"].GetBool():
                        updated_bs_settings.AddString("type", "block_lagrange")
                    else:
                        updated_bs_settings.AddString("type", "block")
                else:
                    updated_bs_settings.AddString("type", "elimination")
            else:
                updated_bs_settings.AddString("type", "block")

            if bs_settings.Has("advanced_settings"):
                updated_bs_settings.AddValue("advanced_settings", bs_settings["advanced_settings"])

            settings.RemoveValue("builder_and_solver_settings")
            settings.AddValue("builder_and_solver_settings", updated_bs_settings)
