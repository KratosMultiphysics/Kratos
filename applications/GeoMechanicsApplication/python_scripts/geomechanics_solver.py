# Importing the Kratos Library
import KratosMultiphysics
from KratosMultiphysics.python_solver import PythonSolver

# Import applications
import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication
import KratosMultiphysics.GeoMechanicsApplication as GeoMechanicsApplication

def CreateSolver(model, custom_settings):
    return GeoMechanicalSolver(model, custom_settings)

class GeoMechanicalSolver(PythonSolver):
    """The base class for geomechanics solvers.

    This class provides shared functions for importing and exporting models,
    adding nodal variables and dofs and solving each solution step.

    Derived classes can override the functions
    """
    def __init__(self, model, custom_settings):

        super().__init__(model, custom_settings)

        self.ValidateSettings()

        model_part_name = self.settings["model_part_name"].GetString()

        if model_part_name == "":
            raise Exception('Please specify a model_part name!')

        # This will be changed once the Model is fully supported!
        if self.model.HasModelPart(model_part_name):
            self.main_model_part = self.model[model_part_name]
            self.solver_imports_model_part = False
        else:
            self.main_model_part = self.model.CreateModelPart(model_part_name)

            domain_size = self.settings["domain_size"].GetInt()
            if domain_size < 0:
                raise Exception('Please specify a "domain_size" >= 0!')

            self.main_model_part.ProcessInfo.SetValue(KratosMultiphysics.DOMAIN_SIZE, domain_size)
            self.solver_imports_model_part = True

        self.min_buffer_size = 2

        KratosMultiphysics.Logger.PrintInfo("::[GeoMechanicalSolver]:: ", "Construction finished")

        # Set if the analysis is restarted
        if self.settings["model_import_settings"]["input_type"].GetString() == "rest":
            KratosMultiphysics.Logger.PrintInfo("geomechanics_solver", "is a restarted model")
            self.main_model_part.ProcessInfo[KratosMultiphysics.IS_RESTARTED] = True
        else:
            KratosMultiphysics.Logger.PrintInfo("geomechanics_solver", "is not a restarted model")
            self.main_model_part.ProcessInfo[KratosMultiphysics.IS_RESTARTED] = False

    @classmethod
    def GetDefaultParameters(cls):
        this_defaults = KratosMultiphysics.Parameters("""{
            "solver_type": "geomechanics_U_Pw_solver",
            "model_part_name": "PorousDomain",
            "domain_size": 2,
            "model_import_settings":{
                "input_type": "mdpa",
                "input_filename": "unknown_name"
            },
            "material_import_settings" :{
                "materials_filename": ""
            },
            "time_stepping": {
                "time_step": 0.1
            },
            "buffer_size": 2,
            "echo_level": 0,
            "rebuild_level": 2,
            "reform_dofs_at_each_step": false,
            "clear_storage": false,
            "compute_reactions": false,
            "move_mesh_flag": false,
            "nodal_smoothing": false,
            "reset_displacements":  false,
            "solution_type": "quasi_static",
            "scheme_type": "Newmark",
            "newmark_beta": 0.25,
            "newmark_gamma": 0.5,
            "newmark_theta": 0.5,
            "rayleigh_m": 0.0,
            "rayleigh_k": 0.0,
            "strategy_type": "newton_raphson",
            "max_piping_iterations": 50,
            "convergence_criterion": "Displacement_criterion",
            "water_pressure_relative_tolerance": 1.0e-4,
            "water_pressure_absolute_tolerance": 1.0e-9,
            "displacement_relative_tolerance": 1.0e-4,
            "displacement_absolute_tolerance": 1.0e-9,
            "residual_relative_tolerance": 1.0e-4,
            "residual_absolute_tolerance": 1.0e-9,
            "desired_iterations"         : 4,
            "max_radius_factor"          : 20.0,
            "min_radius_factor"          : 0.5,
            "max_iterations"             : 15,
            "min_iterations"             : 6,
            "number_cycles"              : 5,
            "increase_factor"            : 2.0,
            "reduction_factor"           : 0.5,
            "calculate_reactions"        : true,
            "max_line_search_iterations" : 5,
            "first_alpha_value"          : 0.5,
            "second_alpha_value"         : 1.0,
            "min_alpha"                  : 0.1,
            "max_alpha"                  : 2.0,
            "line_search_tolerance"      : 0.5,
            "rotation_dofs"              : false,
            "block_builder"              : true,
            "prebuild_dynamics"          : false,
            "search_neighbours_step"     : false,
            "linear_solver_settings":{
                "solver_type": "AMGCL",
                "tolerance": 1.0e-6,
                "max_iteration": 100,
                "scaling": false,
                "verbosity": 0,
                "preconditioner_type": "ILU0Preconditioner",
                "smoother_type": "ilu0",
                "krylov_type": "gmres",
                "coarsening_type": "aggregation"
            },
            "problem_domain_sub_model_part_list": [""],
            "processes_sub_model_part_list": [""],
            "body_domain_sub_model_part_list": [""],
            "loads_sub_model_part_list": [],
            "loads_variable_list": []
        }""")

        this_defaults.AddMissingParameters(super().GetDefaultParameters())
        return this_defaults

    def ValidateSettings(self):
        """This function validates the settings of the solver
        """

        super().ValidateSettings()

        # Checks if scaling is used in combination with rebuild level lower than 2 and prebuild dynamics, if so it
        # throws an error
        if (self.settings.Has("linear_solver_settings") and
            self.settings["linear_solver_settings"].Has("scaling") and
            self.settings["linear_solver_settings"]["scaling"].GetBool()):
            if (self.settings.Has("rebuild_level") and
                self.settings["rebuild_level"].GetInt() < 2):
                raise ValueError("Scaling can only be used if rebuild level is at least equal to 2")
            if (self.settings.Has("prebuild_dynamics") and
                self.settings["prebuild_dynamics"].GetBool()):
                raise ValueError("Scaling can not be used if prebuild dynamics is true")

    def AddVariables(self):
        # this can safely be called also for restarts, it is internally checked if the variables exist already
        # Variables for 1-phase types of calculations:
        # Add displacements.
        self._add_displacement_variables()

        # Add rotational variables
        self._add_rotational_variables()

        # Add dynamic variables
        self._add_dynamic_variables()

        # Variables for 2-phase types of calculations:
        ## Fluid Variables
        self._add_water_variables()

        # Add temperature variables
        self._add_temperature_variables()

        ## smoothing variables
        self._add_smoothing_variables()

        # Add variables that the user defined in the ProjectParameters
        if (self.settings.Has("auxiliary_variables_list")):
            for i in range(self.settings["auxiliary_variables_list"].size()):
                variable_name = self.settings["auxiliary_variables_list"][i].GetString()
                variable = KratosMultiphysics.KratosGlobals.GetVariable(variable_name)
                self.main_model_part.AddNodalSolutionStepVariable(variable)

        KratosMultiphysics.Logger.PrintInfo("::[GeoMechanicalSolver]:: ", "Variables ADDED")

    def AddDofs(self):
        pass

    def GetMinimumBufferSize(self):
        return self.min_buffer_size

    def ImportModelPart(self):
        """This function imports the ModelPart
        """
        if self.solver_imports_model_part:
            self._ImportModelPart(self.main_model_part, self.settings["model_import_settings"])

    def PrepareModelPart(self):
        """This function prepares the ModelPart for being used by the PythonSolver
        """
        # Set ProcessInfo variables
        self.main_model_part.ProcessInfo.SetValue(GeoMechanicsApplication.TIME_UNIT_CONVERTER, 1.0)
        self.main_model_part.ProcessInfo.SetValue(GeoMechanicsApplication.NODAL_SMOOTHING,
                                                  self.settings["nodal_smoothing"].GetBool())

        self.main_model_part.ProcessInfo.SetValue(KratosMultiphysics.STEP, 0)
        if not self.main_model_part.ProcessInfo[KratosMultiphysics.IS_RESTARTED]:
            ## Executes the check and prepare model process (Create computing_model_part and set constitutive law)
            self._ExecuteCheckAndPrepare()
            ## Set buffer size
            self._SetBufferSize()

        if not self.model.HasModelPart(self.settings["model_part_name"].GetString()):
            self.model.AddModelPart(self.main_model_part)

    def KeepAdvancingSolutionLoop(self, end_time):
        current_time_corrected = self.main_model_part.ProcessInfo[KratosMultiphysics.TIME]
        return current_time_corrected < end_time

    def Initialize(self):
        """Perform initialization after adding nodal variables and dofs to the main model part. """
        self.computing_model_part = self.GetComputingModelPart()

        # Fill the previous steps of the buffer with the initial conditions
        self._FillBuffer()

        # Construct the linear solver
        self.linear_solver = self._ConstructLinearSolver()

        # Builder and solver creation
        self.builder_and_solver = self._CreateBuilderAndSolver()

        # Solution scheme creation
        self.scheme = self._ConstructScheme(self.settings["scheme_type"].GetString(),
                                            self.settings["solution_type"].GetString())

        # Get the convergence criterion
        self.convergence_criterion = self._ConstructConvergenceCriterion(self.settings["convergence_criterion"].GetString())

        # Solver creation
        self.solver = self._ConstructSolver(self.builder_and_solver,
                                            self.settings["strategy_type"].GetString())

        # Set echo_level
        self.SetEchoLevel(self.settings["echo_level"].GetInt())

        # Initialize Strategy
        if self.settings["clear_storage"].GetBool():
            self.Clear()

        self.solver.Initialize()

        self.find_neighbour_elements_of_conditions_process = GeoMechanicsApplication.FindNeighbourElementsOfConditionsProcess(self.computing_model_part)
        self.find_neighbour_elements_of_conditions_process.Execute()

        self.deactivate_conditions_on_inactive_elements_process = GeoMechanicsApplication.DeactivateConditionsOnInactiveElements(self.computing_model_part)
        self.deactivate_conditions_on_inactive_elements_process.Execute()

    def InitializeSolutionStep(self):
            self.solver.InitializeSolutionStep()

    def Predict(self):
        self.solver.Predict()

    def SolveSolutionStep(self):
        is_converged = self.solver.SolveSolutionStep()
        return is_converged

    def FinalizeSolutionStep(self):
        self.solver.FinalizeSolutionStep()

    def AdvanceInTime(self, current_time):
        dt = self.ComputeDeltaTime()
        current_time_corrected = self.main_model_part.ProcessInfo[KratosMultiphysics.TIME]
        new_time = current_time_corrected + dt
        self.main_model_part.CloneTimeStep(new_time)
        self.main_model_part.ProcessInfo[KratosMultiphysics.STEP] += 1

        return new_time

    def ComputeDeltaTime(self):
        return self.main_model_part.ProcessInfo[KratosMultiphysics.DELTA_TIME]

    def GetComputingModelPart(self):
        return self.main_model_part.GetSubModelPart(self.computing_model_part_name)

    def ExportModelPart(self):
        name_out_file = self.settings["model_import_settings"]["input_filename"].GetString()+".out"
        file = open(name_out_file + ".mdpa","w")
        file.close()
        KratosMultiphysics.ModelPartIO(name_out_file, KratosMultiphysics.IO.WRITE).WriteModelPart(self.main_model_part)

    def SetEchoLevel(self, level):
        self.solver.SetEchoLevel(level)

    def Clear(self):
        self.solver.Clear()

    def Check(self):
        self.solver.Check()


    #### Specific internal functions ####

    def import_constitutive_laws(self):
        materials_filename = self.settings["material_import_settings"]["materials_filename"].GetString()
        KratosMultiphysics.Logger.PrintInfo("::[GeoMechanicalSolver]:: importing constitutive law", materials_filename)
        if (materials_filename != ""):
            # Add constitutive laws and material properties from json file to model parts.
            material_settings = KratosMultiphysics.Parameters("""{"Parameters": {"materials_filename": ""}} """)
            material_settings["Parameters"]["materials_filename"].SetString(materials_filename)
            KratosMultiphysics.ReadMaterialsUtility(material_settings, self.model)
            materials_imported = True
        else:
            materials_imported = False
        return materials_imported

    #### Private functions ####

    def _add_dynamic_variables(self):
        # For being consistent for Serial and Trilinos
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VELOCITY)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.ACCELERATION)

    def _add_displacement_variables(self):
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)
        self.main_model_part.AddNodalSolutionStepVariable(GeoMechanicsApplication.TOTAL_DISPLACEMENT)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.REACTION)
        self.main_model_part.AddNodalSolutionStepVariable(StructuralMechanicsApplication.POINT_LOAD)
        self.main_model_part.AddNodalSolutionStepVariable(StructuralMechanicsApplication.LINE_LOAD)
        self.main_model_part.AddNodalSolutionStepVariable(StructuralMechanicsApplication.SURFACE_LOAD)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VOLUME_ACCELERATION)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NORMAL_CONTACT_STRESS)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.TANGENTIAL_CONTACT_STRESS)

    def _add_rotational_variables(self):
        if (self.settings.Has("rotation_dofs")):
            if self.settings["rotation_dofs"].GetBool():
                # Add specific variables for the problem (rotation dofs).
                self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.ROTATION)
                self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.REACTION_MOMENT)
                self.main_model_part.AddNodalSolutionStepVariable(StructuralMechanicsApplication.POINT_MOMENT)
                self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.ANGULAR_VELOCITY)
                self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.ANGULAR_ACCELERATION)

    def _add_water_variables(self):
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.WATER_PRESSURE)
        # Add reactions for the water pressure
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.REACTION_WATER_PRESSURE)
        # Add dynamic variables
        self.main_model_part.AddNodalSolutionStepVariable(GeoMechanicsApplication.DT_WATER_PRESSURE)
        # Add variables for the water conditions
        self.main_model_part.AddNodalSolutionStepVariable(GeoMechanicsApplication.NORMAL_FLUID_FLUX)
        # Add variables for the water conditions
        self.main_model_part.AddNodalSolutionStepVariable(GeoMechanicsApplication.HYDRAULIC_DISCHARGE)

    def _add_temperature_variables(self):
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.TEMPERATURE)
        # Add dynamic variables
        self.main_model_part.AddNodalSolutionStepVariable(GeoMechanicsApplication.DT_TEMPERATURE)
        # Add variables for the heat conditions
        self.main_model_part.AddNodalSolutionStepVariable(GeoMechanicsApplication.NORMAL_HEAT_FLUX)
        self.main_model_part.AddNodalSolutionStepVariable(GeoMechanicsApplication.AIR_TEMPERATURE)
        self.main_model_part.AddNodalSolutionStepVariable(GeoMechanicsApplication.SOLAR_RADIATION)
        self.main_model_part.AddNodalSolutionStepVariable(GeoMechanicsApplication.AIR_HUMIDITY)
        self.main_model_part.AddNodalSolutionStepVariable(GeoMechanicsApplication.PRECIPITATION)
        self.main_model_part.AddNodalSolutionStepVariable(GeoMechanicsApplication.WIND_SPEED)

    def _add_smoothing_variables(self):
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NODAL_AREA)
        self.main_model_part.AddNodalSolutionStepVariable(GeoMechanicsApplication.NODAL_CAUCHY_STRESS_TENSOR)
        self.main_model_part.AddNodalSolutionStepVariable(GeoMechanicsApplication.NODAL_DAMAGE_VARIABLE)
        self.main_model_part.AddNodalSolutionStepVariable(GeoMechanicsApplication.NODAL_JOINT_AREA)
        self.main_model_part.AddNodalSolutionStepVariable(GeoMechanicsApplication.NODAL_JOINT_WIDTH)
        self.main_model_part.AddNodalSolutionStepVariable(GeoMechanicsApplication.NODAL_JOINT_DAMAGE)

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

    def _GetLinearSolver(self):
        return self.linear_solver

    def _ExecuteCheckAndPrepare(self):

        self.computing_model_part_name = "porous_computational_model_part"

        # Create list of sub sub model parts (it is a copy of the standard lists with a different name)
        import json

        self.body_domain_sub_sub_model_part_list = []
        for i in range(self.settings["body_domain_sub_model_part_list"].size()):
            self.body_domain_sub_sub_model_part_list.append("sub_"+self.settings["body_domain_sub_model_part_list"][i].GetString())
        self.body_domain_sub_sub_model_part_list = KratosMultiphysics.Parameters(json.dumps(self.body_domain_sub_sub_model_part_list))

        self.loads_sub_sub_model_part_list = []
        for i in range(self.settings["loads_sub_model_part_list"].size()):
            self.loads_sub_sub_model_part_list.append("sub_"+self.settings["loads_sub_model_part_list"][i].GetString())
        self.loads_sub_sub_model_part_list = KratosMultiphysics.Parameters(json.dumps(self.loads_sub_sub_model_part_list))

        # Auxiliary parameters object for the CheckAndPepareModelProcess
        params = KratosMultiphysics.Parameters("{}")
        params.AddEmptyValue("computing_model_part_name").SetString(self.computing_model_part_name)
        params.AddValue("problem_domain_sub_model_part_list",self.settings["problem_domain_sub_model_part_list"])
        params.AddValue("processes_sub_model_part_list",self.settings["processes_sub_model_part_list"])
        params.AddValue("body_domain_sub_model_part_list",self.settings["body_domain_sub_model_part_list"])
        params.AddValue("body_domain_sub_sub_model_part_list",self.body_domain_sub_sub_model_part_list)
        params.AddValue("loads_sub_model_part_list",self.settings["loads_sub_model_part_list"])
        params.AddValue("loads_sub_sub_model_part_list",self.loads_sub_sub_model_part_list)
        # CheckAndPrepareModelProcess creates the porous_computational_model_part
        from KratosMultiphysics.GeoMechanicsApplication import check_and_prepare_model_process_geo
        check_and_prepare_model_process_geo.CheckAndPrepareModelProcess(self.main_model_part, params).Execute()

        # NOTE: We do this here in case the model is empty, so the properties can be assigned
        if not self.model.HasModelPart(self.main_model_part.Name):
            self.model.AddModelPart(self.main_model_part)

        # Import constitutive laws.
        materials_imported = self.import_constitutive_laws()
        if materials_imported:
            KratosMultiphysics.Logger.PrintInfo("::[GeoMechanicalSolver]:: ", "Constitutive law was successfully imported.")
        else:
            raise RuntimeError("::[GeoMechanicalSolver]:: Constitutive law was not imported.")

    def _SetBufferSize(self):
        required_buffer_size = max( self.settings["buffer_size"].GetInt(), self.GetMinimumBufferSize())
        current_buffer_size  = self.main_model_part.GetBufferSize()
        buffer_size          = max(current_buffer_size, required_buffer_size)
        self.main_model_part.SetBufferSize(buffer_size)

    def _FillBuffer(self):
        buffer_size = self.main_model_part.GetBufferSize()
        time        = self.main_model_part.ProcessInfo[KratosMultiphysics.TIME]
        delta_time  = self.main_model_part.ProcessInfo[KratosMultiphysics.DELTA_TIME]
        step        = self.main_model_part.ProcessInfo[KratosMultiphysics.STEP]

        step -= (buffer_size - 1)
        self.main_model_part.ProcessInfo.SetValue(KratosMultiphysics.STEP, step)
        time -= ((buffer_size - 1) * delta_time)
        self.main_model_part.ProcessInfo.SetValue(KratosMultiphysics.TIME, time)
        for _ in range(buffer_size - 1):
            step += 1
            self.main_model_part.ProcessInfo.SetValue(KratosMultiphysics.STEP, step)
            time += delta_time
            self.main_model_part.CloneTimeStep(time)

    def _ConstructLinearSolver(self):
        import KratosMultiphysics.python_linear_solver_factory as linear_solver_factory
        return linear_solver_factory.ConstructSolver(self.settings["linear_solver_settings"])

    def _CreateBuilderAndSolver(self):
        block_builder = self.settings["block_builder"].GetBool()

        # Creating the builder and solver
        if (block_builder):
            builder_and_solver = KratosMultiphysics.ResidualBasedBlockBuilderAndSolver(self.linear_solver)
        else:
            builder_and_solver = KratosMultiphysics.ResidualBasedEliminationBuilderAndSolver(self.linear_solver)

        return builder_and_solver

    def _ConstructSolver(self, builder_and_solver, strategy_type):

        self.main_model_part.ProcessInfo.SetValue(GeoMechanicsApplication.IS_CONVERGED, True)
        self.main_model_part.ProcessInfo.SetValue(KratosMultiphysics.NL_ITERATION_NUMBER, 1)

        max_iters         = self.settings["max_iterations"].GetInt()
        compute_reactions = self.settings["compute_reactions"].GetBool()
        reform_step_dofs  = self.settings["reform_dofs_at_each_step"].GetBool()
        move_mesh_flag    = self.settings["move_mesh_flag"].GetBool()

        if strategy_type.lower() == "newton_raphson":
            self.strategy_params = KratosMultiphysics.Parameters("{}")
            self.strategy_params.AddValue("loads_sub_model_part_list",self.loads_sub_sub_model_part_list)
            self.strategy_params.AddValue("loads_variable_list",self.settings["loads_variable_list"])
            solving_strategy = GeoMechanicsApplication.GeoMechanicsNewtonRaphsonStrategy(self.computing_model_part,
                                                                                         self.scheme,
                                                                                         self.linear_solver,
                                                                                         self.convergence_criterion,
                                                                                         builder_and_solver,
                                                                                         self.strategy_params,
                                                                                         max_iters,
                                                                                         compute_reactions,
                                                                                         reform_step_dofs,
                                                                                         move_mesh_flag)
        elif strategy_type.lower() == "newton_raphson_with_piping":
            self.strategy_params = KratosMultiphysics.Parameters("{}")
            self.strategy_params.AddValue("loads_sub_model_part_list",self.loads_sub_sub_model_part_list)
            self.strategy_params.AddValue("loads_variable_list",self.settings["loads_variable_list"])
            self.strategy_params.AddValue("max_piping_iterations", self.settings["max_piping_iterations"])
            solving_strategy = GeoMechanicsApplication.GeoMechanicsNewtonRaphsonErosionProcessStrategy(self.computing_model_part,
                                                                                                       self.scheme,
                                                                                                       self.linear_solver,
                                                                                                       self.convergence_criterion,
                                                                                                       builder_and_solver,
                                                                                                       self.strategy_params,
                                                                                                       max_iters,
                                                                                                       compute_reactions,
                                                                                                       reform_step_dofs,
                                                                                                       move_mesh_flag)

        elif strategy_type.lower() == "line_search":
            self.strategy_params = KratosMultiphysics.Parameters("{}")
            self.strategy_params.AddValue("max_iteration",              self.settings["max_iterations"])
            self.strategy_params.AddValue("compute_reactions",          self.settings["compute_reactions"])
            self.strategy_params.AddValue("max_line_search_iterations", self.settings["max_line_search_iterations"])
            self.strategy_params.AddValue("first_alpha_value",          self.settings["first_alpha_value"])
            self.strategy_params.AddValue("second_alpha_value",         self.settings["second_alpha_value"])
            self.strategy_params.AddValue("min_alpha",                  self.settings["min_alpha"])
            self.strategy_params.AddValue("max_alpha",                  self.settings["max_alpha"])
            self.strategy_params.AddValue("line_search_tolerance",      self.settings["line_search_tolerance"])
            self.strategy_params.AddValue("move_mesh_flag",             self.settings["move_mesh_flag"])
            self.strategy_params.AddValue("reform_dofs_at_each_step",   self.settings["reform_dofs_at_each_step"])
            self.strategy_params.AddValue("echo_level",                 self.settings["echo_level"])

            solving_strategy = KratosMultiphysics.LineSearchStrategy(self.computing_model_part,
                                                                     self.scheme,
                                                                     self.linear_solver,
                                                                     self.convergence_criterion,
                                                                     self.strategy_params)

        elif strategy_type.lower() == "arc_length":
            # Arc-Length strategy
            self.main_model_part.ProcessInfo.SetValue(KratosGeo.ARC_LENGTH_LAMBDA,        1.0)
            self.main_model_part.ProcessInfo.SetValue(KratosGeo.ARC_LENGTH_RADIUS_FACTOR, 1.0)
            self.strategy_params = KratosMultiphysics.Parameters("{}")
            self.strategy_params.AddValue("desired_iterations",self.settings["desired_iterations"])
            self.strategy_params.AddValue("max_radius_factor",self.settings["max_radius_factor"])
            self.strategy_params.AddValue("min_radius_factor",self.settings["min_radius_factor"])
            self.strategy_params.AddValue("loads_sub_model_part_list",self.loads_sub_sub_model_part_list)
            self.strategy_params.AddValue("loads_variable_list",self.settings["loads_variable_list"])
            solving_strategy = GeoMechanicsApplication.GeoMechanicsRammArcLengthStrategy(self.computing_model_part,
                                                                                         self.scheme,
                                                                                         self.linear_solver,
                                                                                         self.convergence_criterion,
                                                                                         builder_and_solver,
                                                                                         self.strategy_params,
                                                                                         max_iters,
                                                                                         compute_reactions,
                                                                                         reform_step_dofs,
                                                                                         move_mesh_flag)

        elif strategy_type.lower() == "linear":
            solving_strategy = KratosMultiphysics.ResidualBasedLinearStrategy(self.computing_model_part,
                                                                              self.scheme,
                                                                              builder_and_solver,
                                                                              compute_reactions,
                                                                              reform_step_dofs,
                                                                              False,
                                                                              move_mesh_flag)

        else:
            raise RuntimeError(f"Undefined strategy type '{strategy_type}'")

        return solving_strategy

    def _MakeResidualCriterion(self):
        relative_tolerance = self.settings["residual_relative_tolerance"].GetDouble()
        absolute_tolerance = self.settings["residual_absolute_tolerance"].GetDouble()
        residual_criterion = KratosMultiphysics.ResidualCriteria(relative_tolerance, absolute_tolerance)
        residual_criterion.SetEchoLevel(self.settings["echo_level"].GetInt())

        return residual_criterion

