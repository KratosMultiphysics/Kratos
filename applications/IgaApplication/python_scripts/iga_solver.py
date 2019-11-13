from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics

# Import applications
import KratosMultiphysics.IgaApplication as IgaApplication

# Importing the base class
from KratosMultiphysics.python_solver import PythonSolver

# Other imports
from KratosMultiphysics import auxiliary_solver_utilities

def CreateSolver(model, custom_settings):
    return IgaSolver(model, custom_settings)


class IgaSolver(PythonSolver):
    """The base class for iga solvers.

    This class provides functions for importing and exporting models,
    adding nodal variables and dofs and solving each solution step.

    Derived classes must override the function _create_solution_scheme which
    constructs and returns a solution scheme. Depending on the type of
    solver, derived classes may also need to override the following functions:

    _create_solution_scheme
    _create_convergence_criterion
    _create_linear_solver
    _create_builder_and_solver
    _create_solution_strategy

    The iga_solution_strategy, builder_and_solver, etc. should alway be retrieved
    using the getter functions get_solution_strategy, get_builder_and_solver,
    etc. from this base class.

    Only the member variables listed below should be accessed directly.

    Public member variables:
    model -- the model containing the modelpart used to construct the solver.
    settings -- Kratos parameters containing solver settings.
    """
    def __init__(self, model, custom_settings):
        super(IgaSolver, self).__init__(model, custom_settings)

        default_settings = KratosMultiphysics.Parameters("""
        {
            "model_part_name" : "",
            "domain_size" : -1,
            "echo_level": 0,
            "buffer_size": 1,
            "analysis_type": "non_linear",
            "model_import_settings": {
                "input_type": "json",
                "input_filename": "geometry"
            },
            "material_import_settings" :{
                "materials_filename": "materials.json"
            },
            "time_stepping" : { },
            "rotation_dofs": false,
            "reform_dofs_at_each_step": false,
            "line_search": false,
            "compute_reactions": true,
            "block_builder": true,
            "clear_storage": false,
            "move_mesh_flag": true,
            "convergence_criterion": "residual_criterion",
            "displacement_relative_tolerance": 1.0e-4,
            "displacement_absolute_tolerance": 1.0e-9,
            "residual_relative_tolerance": 1.0e-4,
            "residual_absolute_tolerance": 1.0e-9,
            "max_iteration": 10,
            "solver_type": "static",
            "linear_solver_settings":{ },
            "auxiliary_variables_list" : [],
            "auxiliary_dofs_list" : [],
            "auxiliary_reaction_list" : []
        }
        """)

        # Overwrite the default settings with user-provided parameters.
        self.settings.ValidateAndAssignDefaults(default_settings)
        model_part_name = self.settings["model_part_name"].GetString()

        if model_part_name == "":
            raise Exception('Please specify a model_part name!')

        # This will be changed once the Model is fully supported!
        if self.model.HasModelPart(model_part_name):
            self.main_model_part = self.model[model_part_name]
        else:
            self.main_model_part = self.model.CreateModelPart(model_part_name)
            domain_size = self.settings["domain_size"].GetInt()
            if domain_size < 0:
                raise Exception('Please specify a "domain_size" > 0!')
            self.main_model_part.ProcessInfo.SetValue(KratosMultiphysics.DOMAIN_SIZE, domain_size)

        KratosMultiphysics.Logger.PrintInfo("::[IgaSolver]:: ", "Construction finished")

    def AddVariables(self):
        # this can safely be called also for restarts, it is internally checked if the variables exist already
        # Add displacements.
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.REACTION)
        if self.settings["rotation_dofs"].GetBool():
            # Add specific variables for the problem (rotation dofs).
            self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.ROTATION)
            self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.REACTION_MOMENT)
        # Add variables that the user defined in the ProjectParameters
        auxiliary_solver_utilities.AddVariables(self.main_model_part, self.settings["auxiliary_variables_list"])
        KratosMultiphysics.Logger.PrintInfo("::[IgaSolver]:: ", "Variables ADDED")

    def GetMinimumBufferSize(self):
        return max(1, self.settings["buffer_size"].GetDouble())

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
        KratosMultiphysics.Logger.PrintInfo("::[IgaSolver]:: ", "DOF's ADDED")

    def ImportModelPart(self):
        """This function imports the ModelPart"""
        self._set_nurbs_brep_modeler()

        with open(self.settings["model_import_settings"]["physics_filename"].GetString(),'r') as physics_file:
                physics_parameters = KratosMultiphysics.Parameters( physics_file.read())

        physics_file = open(self.settings["model_import_settings"]["physics_filename"].GetString(),'r')
        physics_parameters = KratosMultiphysics.Parameters( physics_file.read())
        self.nurbs_brep_modeler.ImportModelPart(self.main_model_part, physics_parameters)

    def PrepareModelPart(self):
        # Check and prepare computing model part and import constitutive laws.
        self._execute_after_reading()
        self._set_and_fill_buffer()
        KratosMultiphysics.Logger.PrintInfo("::[IgaSolver]::", "ModelPart prepared for Solver.")

    def Initialize(self):
        """Perform initialization after adding nodal variables and dofs to the main model part. """
        KratosMultiphysics.Logger.PrintInfo("::[IgaSolver]:: ", "Initializing ...")
        # The mechanical solution strategy is created here if it does not already exist.
        if self.settings["clear_storage"].GetBool():
            self.Clear()
        iga_solution_strategy = self.get_solution_strategy()
        iga_solution_strategy.SetEchoLevel(self.settings["echo_level"].GetInt())
        iga_solution_strategy.Initialize()

        KratosMultiphysics.Logger.PrintInfo("::[IgaSolver]:: ", "Finished initialization.")

    def InitializeSolutionStep(self):
        if self.settings["clear_storage"].GetBool():
            self.Clear()
            self.Initialize() #required after clearing
        self.get_solution_strategy().InitializeSolutionStep()

    def Predict(self):
        self.get_solution_strategy().Predict()

    def SolveSolutionStep(self):
        is_converged = self.get_solution_strategy().SolveSolutionStep()
        if not is_converged:
            msg  = "Solver did not converge for step " + str(self.main_model_part.ProcessInfo[KratosMultiphysics.STEP]) + "\n"
            msg += "corresponding to time " + str(self.main_model_part.ProcessInfo[KratosMultiphysics.TIME]) + "\n"
            self.print_warning_on_rank_zero("::[IgaSolver]:: ",msg)
        return is_converged

    def FinalizeSolutionStep(self):
        self.get_solution_strategy().FinalizeSolutionStep()

    def AdvanceInTime(self, current_time):
        dt = self._ComputeDeltaTime()
        new_time = current_time + dt
        self.main_model_part.ProcessInfo[KratosMultiphysics.STEP] += 1
        self.main_model_part.CloneTimeStep(new_time)

        return new_time

    def _ComputeDeltaTime(self):
        return self.settings["time_stepping"]["time_step"].GetDouble()

    def GetComputingModelPart(self):
        return self.main_model_part

    def ExportModelPart(self):
        name_out_file = self.settings["model_import_settings"]["input_filename"].GetString()+".out"
        file = open(name_out_file + ".mdpa","w")
        file.close()
        KratosMultiphysics.ModelPartIO(name_out_file, KratosMultiphysics.IO.WRITE).WriteModelPart(self.main_model_part)

    def SetEchoLevel(self, level):
        self.get_solution_strategy().SetEchoLevel(level)

    def Clear(self):
        self.get_solution_strategy().Clear()

    def Check(self):
        self.get_solution_strategy().Check()

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

    def get_solution_strategy(self):
        if not hasattr(self, '_solution_strategy'):
            self._solution_strategy = self._create_solution_strategy()
        return self._solution_strategy

    def import_constitutive_laws(self):
        materials_filename = self.settings["material_import_settings"]["materials_filename"].GetString()
        if (materials_filename != ""):
            # Add constitutive laws and material properties from json file to model parts.
            material_settings = KratosMultiphysics.Parameters("""{"Parameters": {"materials_filename": ""}} """)
            material_settings["Parameters"]["materials_filename"].SetString(materials_filename)
            import read_materials_process # TODO update this, it is deprecated!!!
            read_materials_process.Factory(material_settings,self.model)
            #KratosMultiphysics.ReadMaterialsUtility(material_settings, self.model)
            materials_imported = True
        else:
            materials_imported = False
        return materials_imported

    #### Private functions ####

    def _set_nurbs_brep_modeler(self):
        """Prepare the nurbs brep modeler and read in the necessary data. """
        # This function prepares the nurbs brep modeler and reads in the rough geometry data,
        # which can be used for surface descriptions and integration domains.

        self.nurbs_brep_modeler = IgaApplication.NurbsBrepModeler(self.main_model_part)

        if self.settings["model_import_settings"]["input_type"].GetString() == "json":
            with open(self.settings["model_import_settings"]["input_filename"].GetString() + ".json",'r') as geometry_file:
                geometry_parameters = KratosMultiphysics.Parameters( geometry_file.read())
            self.geometry_reader = IgaApplication.BrepJsonIO()

            self.nurbs_brep_modeler.ImportGeometry(self.geometry_reader, geometry_parameters)

    def _execute_after_reading(self):
        # Import constitutive laws.
        materials_imported = self.import_constitutive_laws()
        if materials_imported:
            self.print_on_rank_zero("::[IgaSolver]:: ", "Constitutive law was successfully imported.")
        else:
            self.print_on_rank_zero("::[IgaSolver]:: ", "Constitutive law was not imported.")

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
        for _ in range(0, buffer_size):
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
        import KratosMultiphysics.StructuralMechanicsApplication.convergence_criteria_factory as convergence_criteria_factory
        convergence_criterion = convergence_criteria_factory.convergence_criterion(self._get_convergence_criterion_settings())
        return convergence_criterion.mechanical_convergence_criterion

    def _create_linear_solver(self):
        linear_solver_configuration = self.settings["linear_solver_settings"]
        if linear_solver_configuration.Has("solver_type"): # user specified a linear solver
            from KratosMultiphysics import python_linear_solver_factory as linear_solver_factory
            return linear_solver_factory.ConstructSolver(linear_solver_configuration)
        else:
            # using a default linear solver (selecting the fastest one available)
            import KratosMultiphysics.kratos_utilities as kratos_utils
            if kratos_utils.CheckIfApplicationsAvailable("EigenSolversApplication"):
                from KratosMultiphysics import EigenSolversApplication
            elif kratos_utils.CheckIfApplicationsAvailable("ExternalSolversApplication"):
                from KratosMultiphysics import ExternalSolversApplication

            linear_solvers_by_speed = [
                "pardiso_lu", # EigenSolversApplication (if compiled with Intel-support)
                "sparse_lu",  # EigenSolversApplication
                "pastix",     # ExternalSolversApplication (if Pastix is included in compilation)
                "super_lu",   # ExternalSolversApplication
                "skyline_lu_factorization" # in Core, always available, but slow
            ]

            for solver_name in linear_solvers_by_speed:
                if KratosMultiphysics.LinearSolverFactory().Has(solver_name):
                    linear_solver_configuration.AddEmptyValue("solver_type").SetString(solver_name)
                    self.print_on_rank_zero('::[MechanicalSolver]:: ',\
                        'Using "' + solver_name + '" as default linear solver')
                    return KratosMultiphysics.LinearSolverFactory().Create(linear_solver_configuration)

        raise Exception("Linear-Solver could not be constructed!")

    def _create_builder_and_solver(self):
        linear_solver = self.get_linear_solver()
        if self.settings["block_builder"].GetBool():
            builder_and_solver = KratosMultiphysics.ResidualBasedBlockBuilderAndSolver(linear_solver)
        else:
            builder_and_solver = KratosMultiphysics.ResidualBasedEliminationBuilderAndSolver(linear_solver)
        return builder_and_solver

    def _create_solution_scheme(self):
        """Create the solution scheme for the structural problem.
        """
        raise Exception("Solution Scheme creation must be implemented in the derived class.")

    def _create_solution_strategy(self):
        analysis_type = self.settings["analysis_type"].GetString()
        if analysis_type == "linear":
            iga_solution_strategy = self._create_linear_strategy()
        elif analysis_type == "non_linear":
            if(self.settings["line_search"].GetBool() == False):
                iga_solution_strategy = self._create_newton_raphson_strategy()
            else:
                iga_solution_strategy = self._create_line_search_strategy()
        else:
            err_msg =  "The requested analysis type \"" + analysis_type + "\" is not available!\n"
            err_msg += "Available options are: \"linear\", \"non_linear\""
            raise Exception(err_msg)
        return iga_solution_strategy

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
