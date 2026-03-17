import KratosMultiphysics

# Import application
import KratosMultiphysics.SPHApplication as SPHApplication

# Import base class 
from KratosMultiphysics.python_solver import PythonSolver

# Other imports
from KratosMultiphysics import python_linear_solver_factory as linear_solver_factory
from KratosMultiphysics.StructuralMechanicsApplication import convergence_criteria_factory
from KratosMultiphysics import auxiliary_solver_utilities
from importlib import import_module

class SPHSolver(PythonSolver):
    """The base class for Smoothed Particle Hydridynamic solvers.   
    
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
    
    The sph_solution_strategy, builder_and_solver, etc. should always be retrieved
    using the getter functions _GetSolutionStrategy, get_builder_and_solver,
    etc. from this base class.  
    
    Only the member variables listed below should be accessed directly. 
    Public member variables:
    model -- the model containing the modelpart used to construct the solver.
    settings -- Kratos parameters containing solver settings.
    """

    def __init__(self, model, custom_settings):
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
        KratosMultiphysics.Logger.PrintInfo("::[SPHSolver]:: ", "Construction finished")

    @classmethod
    def GetDefaultParameters(cls):
        this_defaults = KratosMultiphysics.Parameters("""{
            "solver_type" : "sph_solver",
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
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VOLUME_ACCELERATION)

        # Add variables that the user defined in the ProjectParameters
        auxiliary_solver_utilities.AddVariables(self.main_model_part, self.settings["auxiliary_variables_list"])
        KratosMultiphysics.Logger.PrintInfo("::[SPHSolver]:: ", "Variables ADDED")

    def GetMinimumBufferSize(self):
        return 2

    def AddDofs(self):
        # Append formulation-related DOFs and reactions
        dofs_and_reactions_to_add = []
        dofs_and_reactions_to_add.append(["DISPLACEMENT_X", "REACTION_X"])
        dofs_and_reactions_to_add.append(["DISPLACEMENT_Y", "REACTION_Y"])
        dofs_and_reactions_to_add.append(["DISPLACEMENT_Z", "REACTION_Z"])


        # Append user-defined DOFs and reactions in the ProjectParameters
        auxiliary_solver_utilities.AddAuxiliaryDofsToDofsWithReactionsList(
            self.settings["auxiliary_dofs_list"],
            self.settings["auxiliary_reaction_list"],
            dofs_and_reactions_to_add)

        KratosMultiphysics.VariableUtils.AddDofsList(dofs_and_reactions_to_add, self.main_model_part)
        KratosMultiphysics.Logger.PrintInfo("::[SPHSolver]:: ", "DOF's ADDED")

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

        KratosMultiphysics.Logger.PrintInfo("::[SPHSolver]::", "ModelPart prepared for Solver.")

    def Initialize(self):
        """Perform initialization after adding nodal variables and dofs to the main model part. """
        KratosMultiphysics.Logger.PrintInfo("::[SPHSolver]:: Initializing ...")

        # The mechanical solution strategy is created here if it does not already exist.
        if self.settings["clear_storage"].GetBool():
            self.Clear()
        sph_solution_strategy = self._GetSolutionStrategy()
        sph_solution_strategy.SetEchoLevel(self.settings["echo_level"].GetInt())
        sph_solution_strategy.Initialize()

        # Printing that inialization is finished
        KratosMultiphysics.Logger.PrintInfo("::[SPHSolver]:: Finished initialization.")

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
            KratosMultiphysics.Logger.PrintWarning("::[SPHSolver]:: ",msg)
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
            raise Exception("::[SPHSolver]:: Time stepping not defined!")

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
        if not hasattr(self, '_builder_and_solver'):
            self._builder_and_solver = self._CreateBuilderAndSolver()
        return self._builder_and_solver

    def _GetSolutionStrategy(self):
        if not hasattr(self, '_sph_solution_strategy'):
            self._sph_solution_strategy = self._CreateSolutionStrategy()
        return self._sph_solution_strategy

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
            KratosMultiphysics.Logger.PrintInfo("::[SPHSolver]:: ", "Constitutive law was successfully imported.")
        else:
            KratosMultiphysics.Logger.PrintInfo("::[SPHSolver]:: ", "Constitutive law was not imported.")

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

    def _add_dynamic_dofs(self):
        # For being consistent for Serial and Trilinos
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.VELOCITY_X,self.main_model_part)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.VELOCITY_Y,self.main_model_part)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.VELOCITY_Z,self.main_model_part)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.ACCELERATION_X,self.main_model_part)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.ACCELERATION_Y,self.main_model_part)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.ACCELERATION_Z,self.main_model_part)

    def _get_convergence_criterion_settings(self):
        # Create an auxiliary Kratos parameters object to store the convergence settings.
        conv_params = KratosMultiphysics.Parameters("{}")
        conv_params.AddValue("convergence_criterion",self.settings["convergence_criterion"])
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
            KratosMultiphysics.Logger.PrintInfo('::[SPHSolver]:: No linear solver was specified, using fastest available solver')
            return linear_solver_factory.CreateFastestAvailableDirectLinearSolver()

    def _CreateBuilderAndSolver(self) -> KratosMultiphysics.BuilderAndSolver:
        if not self.settings["builder_and_solver_settings"].Has("type"):
            return self.__DeprecatedBuilderAndSolverFactory()

        linear_solver: KratosMultiphysics.LinearSolver = self._GetLinearSolver()
        builder_and_solver_name: str = self.settings["builder_and_solver_settings"]["type"].GetString()

        if builder_and_solver_name == "block":
            bs_params = self.settings["builder_and_solver_settings"]["advanced_settings"]
            return KratosMultiphysics.ResidualBasedBlockBuilderAndSolver(linear_solver, bs_params)
        else:
            message: str = f"Invalid type for builder and solver '{builder_and_solver_name}'. Options are:\n"
            message += "'block'\n"
            raise ValueError(message)

    def _CreateScheme(self):
        """Create the solution scheme for the structural problem.
        """
        raise Exception("Solution Scheme creation must be implemented in the derived class.")

    def _CreateSolutionStrategy(self):
        analysis_type = self.settings["analysis_type"].GetString()
        if analysis_type == "linear":
            sph_solution_strategy = self._create_linear_strategy()
        elif analysis_type == "non_linear":
            # Create strategy
            if self.settings["solving_strategy_settings"]["type"].GetString() == "newton_raphson":
                sph_solution_strategy = self._create_newton_raphson_strategy()
        else:
            err_msg =  "The requested analysis type \"" + analysis_type + "\" is not available!\n"
            err_msg += "Available options are: \"linear\", \"non_linear\""
            raise Exception(err_msg)
        return sph_solution_strategy

    def _create_linear_strategy(self):
        computing_model_part = self.GetComputingModelPart()
        sph_scheme = self._GetScheme()
        builder_and_solver = self._GetBuilderAndSolver()
        return KratosMultiphysics.ResidualBasedLinearStrategy(computing_model_part,
                                                              sph_scheme,
                                                              builder_and_solver,
                                                              self.settings["compute_reactions"].GetBool(),
                                                              self.settings["reform_dofs_at_each_step"].GetBool(),
                                                              False,
                                                              self.settings["move_mesh_flag"].GetBool())

    def _create_newton_raphson_strategy(self):
        computing_model_part = self.GetComputingModelPart()
        sph_scheme = self._GetScheme()
        sph_convergence_criterion = self._GetConvergenceCriterion()
        builder_and_solver = self._GetBuilderAndSolver()
        strategy = KratosMultiphysics.ResidualBasedNewtonRaphsonStrategy(computing_model_part,
                                                                     sph_scheme,
                                                                     sph_convergence_criterion,
                                                                     builder_and_solver,
                                                                     self.settings["max_iteration"].GetInt(),
                                                                     self.settings["compute_reactions"].GetBool(),
                                                                     self.settings["reform_dofs_at_each_step"].GetBool(),
                                                                     self.settings["move_mesh_flag"].GetBool())
        strategy.SetUseOldStiffnessInFirstIterationFlag(self.settings["use_old_stiffness_in_first_iteration"].GetBool())
        return strategy
