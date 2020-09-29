from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics
from KratosMultiphysics.python_solver import PythonSolver

# Check that applications were imported in the main script
# DEPRECATION: "CheckRegisteredApplications" is not needed any more and can be safely removed
#KratosMultiphysics.CheckRegisteredApplications("StructuralMechanicsApplication")
#KratosMultiphysics.CheckRegisteredApplications("GeoMechanicsApplication")

# Import applications
import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication
import KratosMultiphysics.GeoMechanicsApplication as GeoMechanicsApplication


def CreateSolver(model, custom_settings):
    return GeoMechanicalSolver(model, custom_settings)

class GeoMechanicalSolver(PythonSolver):
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
        self._validate_settings_in_baseclass=False # To be removed eventually
        super(GeoMechanicalSolver, self).__init__(model, custom_settings)

        # default_settings = KratosMultiphysics.Parameters("""
        # {
        #     "model_part_name" : "",
        #     "domain_size" : -1,
        #     "echo_level": 0,
        #     "buffer_size": 2,
        #     "analysis_type": "non_linear",
        #     "model_import_settings": {
        #         "input_type": "mdpa",
        #         "input_filename": "unknown_name"
        #     },
        #     "computing_model_part_name" : "computing_domain",
        #     "material_import_settings" :{
        #         "materials_filename": ""
        #     },
        #     "time_stepping" : { },
        #     "rotation_dofs": false,
        #     "reform_dofs_at_each_step": false,
        #     "line_search": false,
        #     "compute_reactions": true,
        #     "block_builder": true,
        #     "clear_storage": false,
        #     "move_mesh_flag": true,
        #     "multi_point_constraints_used": false,
        #     "convergence_criterion": "residual_criterion",
        #     "displacement_relative_tolerance": 1.0e-4,
        #     "displacement_absolute_tolerance": 1.0e-9,
        #     "residual_relative_tolerance": 1.0e-4,
        #     "residual_absolute_tolerance": 1.0e-9,
        #     "max_iteration": 10,
        #     "linear_solver_settings":{
        #         "solver_type": "SuperLUSolver",
        #         "max_iteration": 500,
        #         "tolerance": 1e-9,
        #         "scaling": false,
        #         "verbosity": 1
        #     },
        #     "problem_domain_sub_model_part_list": ["solid"],
        #     "processes_sub_model_part_list": [""],
        #     "auxiliary_variables_list" : [],
        #     "auxiliary_dofs_list" : [],
        #     "auxiliary_reaction_list" : []
        # }
        # """)

        # # temporary warnings, to be removed
        # if custom_settings.Has("bodies_list"):
        #     custom_settings.RemoveValue("bodies_list")
        #     warning = '\n::[GeoMechanicalSolver]:: W-A-R-N-I-N-G: You have specified "bodies_list", '
        #     warning += 'which is deprecated and will be removed soon. \nPlease remove it from the "solver settings"!\n'
        #     self.print_warning_on_rank_zero("Bodies list", warning)
        # if custom_settings.Has("solver_type"):
        #     custom_settings.RemoveValue("solver_type")
        #     warning = '\n::[GeoMechanicalSolver]:: W-A-R-N-I-N-G: You have specified "solver_type", '
        #     warning += 'which is only needed if you use the "python_solvers_wrapper_structural". \nPlease remove it '
        #     warning += 'from the "solver settings" if you dont use this wrapper, this check will be removed soon!\n'
        #     self.print_warning_on_rank_zero("Solver type", warning)
        # if custom_settings.Has("time_integration_method"):
        #     custom_settings.RemoveValue("time_integration_method")
        #     warning = '\n::[GeoMechanicalSolver]:: W-A-R-N-I-N-G: You have specified "time_integration_method", '
        #     warning += 'which is only needed if you use the "python_solvers_wrapper_structural". \nPlease remove it '
        #     warning += 'from the "solver settings" if you dont use this wrapper, this check will be removed soon!\n'
        #     self.print_warning_on_rank_zero("Time integration method", warning)

        # # Overwrite the default settings with user-provided parameters.
        # self.settings.ValidateAndAssignDefaults(default_settings)
        model_part_name = self.settings["model_part_name"].GetString()

        if model_part_name == "":
            raise Exception('Please specify a model_part name!')

        # This will be changed once the Model is fully supported!
        if self.model.HasModelPart(model_part_name):
            self.main_model_part = self.model[model_part_name]
            self.solver_imports_model_part = False
        else:
            #self.main_model_part = KratosMultiphysics.ModelPart(model_part_name) # Model.CreateodelPart()
            self.main_model_part = self.model.CreateModelPart(model_part_name)
            domain_size = self.settings["domain_size"].GetInt()
            if domain_size < 0:
                raise Exception('Please specify a "domain_size" >= 0!')
            self.main_model_part.ProcessInfo.SetValue(KratosMultiphysics.DOMAIN_SIZE, domain_size)
            self.solver_imports_model_part = True

        KratosMultiphysics.Logger.PrintInfo("::[GeoMechanicalSolver]:: ", "Construction finished")

        # Set if the analysis is restarted
        if self.settings["model_import_settings"]["input_type"].GetString() == "rest":
            KratosMultiphysics.Logger.PrintInfo("geomechanics_solver", "is a restarted model")
            self.main_model_part.ProcessInfo[KratosMultiphysics.IS_RESTARTED] = True
        else:
            KratosMultiphysics.Logger.PrintInfo("geomechanics_solver", "is not a restarted model")
            self.main_model_part.ProcessInfo[KratosMultiphysics.IS_RESTARTED] = False

    def AddVariables(self):
        # this can safely be called also for restarts, it is internally checked if the variables exist already
        # Variables for 1-phase types of calculations:
        # Add displacements.
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)
        self.main_model_part.AddNodalSolutionStepVariable(GeoMechanicsApplication.TOTAL_DISPLACEMENT)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.REACTION)
        # Add specific variables for the problem conditions.
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.POSITIVE_FACE_PRESSURE)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NEGATIVE_FACE_PRESSURE)
        self.main_model_part.AddNodalSolutionStepVariable(StructuralMechanicsApplication.POINT_LOAD)
        self.main_model_part.AddNodalSolutionStepVariable(StructuralMechanicsApplication.LINE_LOAD)
        self.main_model_part.AddNodalSolutionStepVariable(StructuralMechanicsApplication.SURFACE_LOAD)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VOLUME_ACCELERATION)
        if (self.settings.Has("rotation_dofs")):
            if self.settings["rotation_dofs"].GetBool():
                # Add specific variables for the problem (rotation dofs).
                self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.ROTATION)
                self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.REACTION_MOMENT)
                self.main_model_part.AddNodalSolutionStepVariable(StructuralMechanicsApplication.POINT_MOMENT)
                self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.ANGULAR_VELOCITY)
                self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.ANGULAR_ACCELERATION)

        # Add variables that the user defined in the ProjectParameters
        if (self.settings.Has("auxiliary_variables_list")):
            for i in range(self.settings["auxiliary_variables_list"].size()):
                variable_name = self.settings["auxiliary_variables_list"][i].GetString()
                variable = KratosMultiphysics.KratosGlobals.GetVariable(variable_name)
                self.main_model_part.AddNodalSolutionStepVariable(variable)

        # Variables for 2-phase types of calculations:
        # Add dynamic variables
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VELOCITY)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.ACCELERATION)
        # Add variables for the solid conditions
        self.main_model_part.AddNodalSolutionStepVariable(StructuralMechanicsApplication.POINT_LOAD)
        self.main_model_part.AddNodalSolutionStepVariable(StructuralMechanicsApplication.LINE_LOAD)
        self.main_model_part.AddNodalSolutionStepVariable(StructuralMechanicsApplication.SURFACE_LOAD)
        #self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.FACE_LOAD)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NORMAL_CONTACT_STRESS)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.TANGENTIAL_CONTACT_STRESS)
        ## Fluid Variables
        # Add water pressure
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.WATER_PRESSURE)
        # Add reactions for the water pressure
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.REACTION_WATER_PRESSURE)
        # Add dynamic variables
        self.main_model_part.AddNodalSolutionStepVariable(GeoMechanicsApplication.DT_WATER_PRESSURE)
        # Add variables for the water conditions
        self.main_model_part.AddNodalSolutionStepVariable(GeoMechanicsApplication.NORMAL_FLUID_FLUX)
        ## Other variables
        if (self.settings.Has("nodal_smoothing")):
            if(self.settings["nodal_smoothing"].GetBool() == True):
                self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NODAL_AREA)
                self.main_model_part.AddNodalSolutionStepVariable(GeoMechanicsApplication.NODAL_CAUCHY_STRESS_TENSOR)
                self.main_model_part.AddNodalSolutionStepVariable(GeoMechanicsApplication.NODAL_DAMAGE_VARIABLE)
                self.main_model_part.AddNodalSolutionStepVariable(GeoMechanicsApplication.NODAL_JOINT_AREA)
                self.main_model_part.AddNodalSolutionStepVariable(GeoMechanicsApplication.NODAL_JOINT_WIDTH)
                self.main_model_part.AddNodalSolutionStepVariable(GeoMechanicsApplication.NODAL_JOINT_DAMAGE)

        KratosMultiphysics.Logger.PrintInfo("::[GeoMechanicalSolver]:: ", "Variables ADDED")

    def GetMinimumBufferSize(self):
        return 2

    def AddDofs(self):
        pass


    def ImportModelPart(self):
        """This function imports the ModelPart
        """
        if self.solver_imports_model_part:
            self._ImportModelPart(self.main_model_part, self.settings["model_import_settings"])

    def PrepareModelPart(self):
        """This function prepares the ModelPart for being used by the PythonSolver
        """
        pass

    def KeepAdvancingSolutionLoop(self):
        pass

    # def PrepareModelPart(self):
    #     if not self.is_restarted():
    #         # Check and prepare computing model part and import constitutive laws.
    #         self._execute_after_reading()
    #         self._set_and_fill_buffer()

    #     # This will be removed once the Model is fully supported! => It wont e necessary anymore
    #     if not self.model.HasModelPart(self.main_model_part.Name):
    #         self.model.AddModelPart(self.main_model_part)

    #     KratosMultiphysics.Logger.PrintInfo("::[GeoMechanicalSolver]::", "ModelPart prepared for Solver.")

    def Initialize(self):
        """Perform initialization after adding nodal variables and dofs to the main model part. """
        KratosMultiphysics.Logger.PrintInfo("::[GeoMechanicalSolver]:: ", "Initializing ...")
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
        self.Check()
        KratosMultiphysics.Logger.PrintInfo("::[GeoMechanicalSolver]:: ", "Finished initialization.")

    # def Solve(self):
    #     KratosMultiphysics.Logger.PrintWarning("::[GeoMechanicalSolver]:: ", "Solve(self)")
    #     if self.settings["clear_storage"].GetBool():
    #         self.Clear()
    #     mechanical_solution_strategy = self.get_mechanical_solution_strategy()
    #     mechanical_solution_strategy.Solve()

    def InitializeSolutionStep(self):
        self.get_mechanical_solution_strategy().InitializeSolutionStep()

    def Predict(self):
        self.get_mechanical_solution_strategy().Predict()

    def SolveSolutionStep(self):
        is_converged = self.get_mechanical_solution_strategy().SolveSolutionStep()
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
        if not self.main_model_part.HasSubModelPart(self.settings["computing_model_part_name"].GetString()):
            raise Exception("The ComputingModelPart was not created yet!")
        return self.main_model_part.GetSubModelPart(self.settings["computing_model_part_name"].GetString())

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
        if not hasattr(self, '_builder_and_solver'):
            self._builder_and_solver = self._create_builder_and_solver()
        return self._builder_and_solver

    def get_mechanical_solution_strategy(self):
        KratosMultiphysics.Logger.PrintInfo("[GeoMechanicalSolver]", "get_mechanical_solution_strategy")
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
        import convergence_criteria_factory
        convergence_criterion = convergence_criteria_factory.convergence_criterion(self._get_convergence_criterion_settings())
        return convergence_criterion.mechanical_convergence_criterion

    def _create_linear_solver(self):
        import python_linear_solver_factory
        linear_solver = python_linear_solver_factory.ConstructSolver(self.settings["linear_solver_settings"])
        return linear_solver

    def _create_builder_and_solver(self):
        linear_solver = self.get_linear_solver()
        if self.settings["block_builder"].GetBool():
            if self.settings["multi_point_constraints_used"].GetBool():
                builder_and_solver = KratosMultiphysics.StructuralMechanicsApplication.ResidualBasedBlockBuilderAndSolverWithMpc(linear_solver)
            else:
                builder_and_solver = KratosMultiphysics.ResidualBasedBlockBuilderAndSolver(linear_solver)
        else:
            if self.settings["multi_point_constraints_used"].GetBool():
                raise Exception("To use MPCs you also have to set \"block_builder\" to \"true\"")
            builder_and_solver = KratosMultiphysics.ResidualBasedEliminationBuilderAndSolver(linear_solver)
        return builder_and_solver

    def _create_solution_scheme(self):
        """Create the solution scheme for the structural problem.
        """
        raise Exception("Solution Scheme creation must be implemented in the derived class.")

    def _create_mechanical_solution_strategy(self):
        analysis_type = self.settings["analysis_type"].GetString()
        KratosMultiphysics.Logger.PrintInfo("[GeoMechanicalSolver]: analysis_type", analysis_type)

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
