from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics

# Import applications and dependencies
import KratosMultiphysics.ParticleMechanicsApplication as KratosParticle

# Importing the base class
from KratosMultiphysics.python_solver import PythonSolver

def CreateSolver(model, custom_settings):
    return MPMSolver(model, custom_settings)

class MPMSolver(PythonSolver):

    ### Solver constructor
    def __init__(self, model, custom_settings):
        self._validate_settings_in_baseclass=True # To be removed eventually
        super(MPMSolver, self).__init__(model, custom_settings)

        # Add model part containers
        self._add_model_part_containers()

        # Default settings
        self.min_buffer_size = 2

        KratosMultiphysics.Logger.PrintInfo("::[MPMSolver]:: ", "Solver is constructed correctly.")


    @classmethod
    def GetDefaultSettings(cls):
        this_defaults = KratosMultiphysics.Parameters("""
        {
            "model_part_name" : "MPM_Material",
            "domain_size"     : -1,
            "echo_level"      : 0,
            "time_stepping"   : { },
            "analysis_type"   : "non_linear",
            "grid_model_import_settings" : {
                "input_type"     : "mdpa",
                "input_filename" : "unknown_name_Grid"
            },
            "model_import_settings" : {
                "input_type"        : "mdpa",
                "input_filename"    : "unknown_name_Body"
            },
            "material_import_settings" : {
                "materials_filename" : ""
            },
            "compute_reactions"                  : false,
            "convergence_criterion"              : "Residual_criteria",
            "displacement_relative_tolerance"    : 1.0E-4,
            "displacement_absolute_tolerance"    : 1.0E-9,
            "residual_relative_tolerance"        : 1.0E-4,
            "residual_absolute_tolerance"        : 1.0E-9,
            "max_iteration"                      : 20,
            "pressure_dofs"                      : false,
            "axis_symmetric_flag"                : false,
            "block_builder"                      : true,
            "move_mesh_flag"                     : false,
            "problem_domain_sub_model_part_list" : [],
            "processes_sub_model_part_list"      : [],
            "auxiliary_variables_list"           : [],
            "auxiliary_dofs_list"                : [],
            "auxiliary_reaction_list"            : [],
            "element_search_settings"            : {
                "search_algorithm_type"          : "bin_based",
                "max_number_of_results"          : 1000,
                "searching_tolerance"            : 1.0E-5
            },
            "linear_solver_settings"             : {
                "solver_type" : "amgcl",
                "smoother_type":"damped_jacobi",
                "krylov_type": "cg",
                "coarsening_type": "aggregation",
                "max_iteration": 200,
                "provide_coordinates": false,
                "gmres_krylov_space_dimension": 100,
                "verbosity" : 0,
                "tolerance": 1e-7,
                "scaling": false,
                "block_size": 3,
                "use_block_matrices_if_possible" : true,
                "coarse_enough" : 50
            }
        }""")
        this_defaults.AddMissingParameters(super(MPMSolver, cls).GetDefaultSettings())
        return this_defaults

    ### Solver public functions

    def AddVariables(self):
        # Add variables to background grid model part
        self._add_variables_to_model_part(self.grid_model_part)

        KratosMultiphysics.Logger.PrintInfo("::[MPMSolver]:: ", "Variables are added.")

    def ImportModelPart(self):
        # Read model part
        self._model_part_reading()

        KratosMultiphysics.Logger.PrintInfo("::[MPMSolver]:: ","Models are imported.")

    def PrepareModelPart(self):
        # Set buffer size
        self._set_buffer_size()

        # Executes the check and prepare model process
        self._execute_check_and_prepare()

        KratosMultiphysics.Logger.PrintInfo("::[MPMSolver]:: ", "ModelPart prepared for Solver.")

    def GetComputingModelPart(self):
        if not self.model.HasModelPart(self.settings["model_part_name"].GetString()):
            raise Exception("The ComputingModelPart was not created yet!")
        return self.model.GetModelPart(self.settings["model_part_name"].GetString())

    def GetGridModelPart(self):
        if not self.model.HasModelPart("Background_Grid"):
            raise Exception("The GridModelPart was not created yet!")
        return self.model.GetModelPart("Background_Grid")

    def AddDofs(self):
        # Add dofs to background grid model part
        self._add_dofs_to_model_part(self.grid_model_part)

        KratosMultiphysics.Logger.PrintInfo("::[MPMSolver]:: ","DOFs are added.")

    def Initialize(self):
        # The mechanical solution strategy is created here if it does not already exist.
        particle_solution_strategy = self.get_solution_strategy()
        particle_solution_strategy.SetEchoLevel(self.settings["echo_level"].GetInt())

        # Generate material points
        self.generate_material_point()

        KratosMultiphysics.Logger.PrintInfo("::[MPMSolver]:: ","Solver is initialized correctly.")

    def AdvanceInTime(self, current_time):
        dt = self.ComputeDeltaTime()
        new_time = current_time + dt
        self.grid_model_part.ProcessInfo[KratosMultiphysics.STEP] += 1
        self.grid_model_part.CloneTimeStep(new_time)

        return new_time

    def ComputeDeltaTime(self):
        return self.settings["time_stepping"]["time_step"].GetDouble()

    def InitializeSolutionStep(self):
        self.search_element()
        self.get_solution_strategy().Initialize()
        self.get_solution_strategy().InitializeSolutionStep()

    def Predict(self):
        self.get_solution_strategy().Predict()

    def SolveSolutionStep(self):
        is_converged = self.get_solution_strategy().SolveSolutionStep()
        return is_converged

    def FinalizeSolutionStep(self):
        self.get_solution_strategy().FinalizeSolutionStep()

        self.get_solution_strategy().Clear()

    def Check(self):
        self.get_solution_strategy().Check()

    def Clear(self):
        self.get_solution_strategy().Clear()

    ### Solver special private functions

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

    def generate_material_point(self):
        pressure_dofs          = self.settings["pressure_dofs"].GetBool()
        axis_symmetric_flag    = self.settings["axis_symmetric_flag"].GetBool()

        # Assigning extra information to the main model part
        self.material_model_part.SetNodes(self.grid_model_part.GetNodes())
        self.material_model_part.ProcessInfo = self.grid_model_part.ProcessInfo
        self.material_model_part.SetBufferSize(self.grid_model_part.GetBufferSize())

        # Generate MP Element and Condition
        KratosParticle.GenerateMaterialPointElement(self.grid_model_part, self.initial_material_model_part, self.material_model_part, axis_symmetric_flag, pressure_dofs)
        KratosParticle.GenerateMaterialPointCondition(self.grid_model_part, self.initial_material_model_part, self.material_model_part)

    def search_element(self):
        searching_alg_type = self.settings["element_search_settings"]["search_algorithm_type"].GetString()
        max_number_of_search_results = self.settings["element_search_settings"]["max_number_of_results"].GetInt()
        searching_tolerance          = self.settings["element_search_settings"]["searching_tolerance"].GetDouble()
        if (searching_alg_type == "bin_based"):
            KratosParticle.SearchElement(self.grid_model_part, self.material_model_part, max_number_of_search_results, searching_tolerance)
        else:
            err_msg  = "The requested searching algorithm \"" + searching_alg_type
            err_msg += "\" is not available for ParticleMechanicsApplication!\n"
            err_msg += "Available options are: \"bin_based\""
            raise Exception(err_msg)


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

    ### Solver private functions

    def _add_model_part_containers(self):

        domain_size = self._get_domain_size()
        if domain_size not in [2,3]:
            err_msg  = "The input \"domain_size\" is wrong!"
            err_msg += "Available options are: \"2\" or \"3\""
            raise Exception(err_msg)

        ### In MPM three model parts are needed
        ## Material model part definition
        material_model_part_name = self.settings["model_part_name"].GetString()
        if not self.model.HasModelPart(material_model_part_name):
            self.material_model_part = self.model.CreateModelPart(material_model_part_name) # Equivalent to model_part3 in the old format
            self.material_model_part.ProcessInfo.SetValue(KratosMultiphysics.DOMAIN_SIZE, domain_size)

        ## Initial material model part definition
        initial_material_model_part_name = "Initial_" + material_model_part_name
        if not self.model.HasModelPart(initial_material_model_part_name):
            self.initial_material_model_part = self.model.CreateModelPart(initial_material_model_part_name) #Equivalent to model_part2 in the old format
            self.initial_material_model_part.ProcessInfo.SetValue(KratosMultiphysics.DOMAIN_SIZE, domain_size)

        ## Grid model part definition
        if not self.model.HasModelPart("Background_Grid"):
            self.grid_model_part = self.model.CreateModelPart("Background_Grid") #Equivalent to model_part1 in the old format
            self.grid_model_part.ProcessInfo.SetValue(KratosMultiphysics.DOMAIN_SIZE, domain_size)

    def _add_variables_to_model_part(self, model_part):
        # Add displacements and reaction
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.REACTION)

        # Add specific variables for the problem conditions
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.PRESSURE)
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.POSITIVE_FACE_PRESSURE)
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VOLUME_ACCELERATION)

        # MPM specific nodal variables
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NODAL_MASS)
        model_part.AddNodalSolutionStepVariable(KratosParticle.NODAL_MOMENTUM)
        model_part.AddNodalSolutionStepVariable(KratosParticle.NODAL_INERTIA)

        # Add variables that the user defined in the ProjectParameters
        for i in range(self.settings["auxiliary_variables_list"].size()):
            variable_name = self.settings["auxiliary_variables_list"][i].GetString()
            variable = KratosMultiphysics.KratosGlobals.GetVariable(variable_name)
            model_part.AddNodalSolutionStepVariable(variable)

        # Add variables for specific cases
        if self.settings["pressure_dofs"].GetBool():
            # add specific variables for the problem (pressure dofs)
            model_part.AddNodalSolutionStepVariable(KratosParticle.PRESSURE_REACTION)
            model_part.AddNodalSolutionStepVariable(KratosParticle.NODAL_MPRESSURE)

    def _add_dynamic_variables(self, model_part):
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VELOCITY)
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.ACCELERATION)

    def _model_part_reading(self):
        # reading the model part of the background grid
        if(self.settings["grid_model_import_settings"]["input_type"].GetString() == "mdpa"):
            self._ImportModelPart(self.grid_model_part, self.settings["grid_model_import_settings"])
        else:
            raise Exception("Other input options are not implemented yet.")

        # reading the model part of the material point
        if(self.settings["model_import_settings"]["input_type"].GetString() == "mdpa"):
            self._ImportModelPart(self.initial_material_model_part, self.settings["model_import_settings"])
        else:
            raise Exception("Other input options are not implemented yet.")

    def _execute_check_and_prepare(self):
        # Specific active node and element check for particle MPM solver
        for node in self.grid_model_part.Nodes:
            if (node.Is(KratosMultiphysics.ACTIVE)):
                KratosMultiphysics.Logger.PrintInfo("::[MPMSolver]:: ","WARNING: This grid node have been set active: ", node.Id)

        # Setting active initial elements
        KratosMultiphysics.VariableUtils().SetFlag(KratosMultiphysics.ACTIVE, True, self.initial_material_model_part.Elements)

        # Read material property
        materials_imported = self.import_constitutive_laws()
        if materials_imported:
            KratosMultiphysics.Logger.PrintInfo("::[MPMSolver]:: ","Constitutive law was successfully imported.")
        else:
            KratosMultiphysics.Logger.PrintWarning("::[MPMSolver]:: ","Constitutive law was not imported.")

        # Clone property of model_part2 to model_part3
        self.material_model_part.Properties = self.initial_material_model_part.Properties

    def _add_dofs_to_model_part(self, model_part):
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_X, KratosMultiphysics.REACTION_X, model_part)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_Y, KratosMultiphysics.REACTION_Y, model_part)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_Z, KratosMultiphysics.REACTION_Z, model_part)

        if self.settings["pressure_dofs"].GetBool():
            KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.PRESSURE, KratosParticle.PRESSURE_REACTION, model_part)

        # Add dofs that the user defined in the ProjectParameters
        if (self.settings["auxiliary_dofs_list"].size() != self.settings["auxiliary_reaction_list"].size()):
                raise Exception("DoFs list and reaction list should be the same long")
        for i in range(self.settings["auxiliary_dofs_list"].size()):
            dof_variable_name = self.settings["auxiliary_dofs_list"][i].GetString()
            reaction_variable_name = self.settings["auxiliary_reaction_list"][i].GetString()
            if (KratosMultiphysics.KratosGlobals.HasDoubleVariable(dof_variable_name)): # Double variable
                dof_variable = KratosMultiphysics.KratosGlobals.GetVariable(dof_variable_name)
                reaction_variable = KratosMultiphysics.KratosGlobals.GetVariable(reaction_variable_name)
                KratosMultiphysics.VariableUtils().AddDof(dof_variable, reaction_variable, model_part)
            elif (KratosMultiphysics.KratosGlobals.HasArrayVariable(dof_variable_name)): # Components variable
                dof_variable_x = KratosMultiphysics.KratosGlobals.GetVariable(dof_variable_name + "_X")
                reaction_variable_x = KratosMultiphysics.KratosGlobals.GetVariable(reaction_variable_name + "_X")
                KratosMultiphysics.VariableUtils().AddDof(dof_variable_x, reaction_variable_x, model_part)
                dof_variable_y = KratosMultiphysics.KratosGlobals.GetVariable(dof_variable_name + "_Y")
                reaction_variable_y = KratosMultiphysics.KratosGlobals.GetVariable(reaction_variable_name + "_Y")
                KratosMultiphysics.VariableUtils().AddDof(dof_variable_y, reaction_variable_y, model_part)
                dof_variable_z = KratosMultiphysics.KratosGlobals.GetVariable(dof_variable_name + "_Z")
                reaction_variable_z = KratosMultiphysics.KratosGlobals.GetVariable(reaction_variable_name + "_Z")
                KratosMultiphysics.VariableUtils().AddDof(dof_variable_z, reaction_variable_z, model_part)
            else:
                KratosMultiphysics.Logger.PrintWarning("auxiliary_reaction_list list", "The variable " + dof_variable_name + "is not a compatible type")

    def _get_domain_size(self):
        if not hasattr(self, '_domain_size'):
            self._domain_size = self.settings["domain_size"].GetInt()
        return self._domain_size

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

    def _create_convergence_criterion(self):
        convergence_criterion_parameters = self._get_convergence_criterion_settings()
        if (convergence_criterion_parameters["convergence_criterion"].GetString() == "residual_criterion"):
            R_RT = convergence_criterion_parameters["residual_relative_tolerance"].GetDouble()
            R_AT = convergence_criterion_parameters["residual_absolute_tolerance"].GetDouble()
            convergence_criterion = KratosMultiphysics.ResidualCriteria(R_RT, R_AT)
            convergence_criterion.SetEchoLevel(convergence_criterion_parameters["echo_level"].GetInt())
        else:
            err_msg  = "The requested convergence criteria \"" + convergence_criterion_parameters["convergence_criterion"].GetString()
            err_msg += "\" is not supported for ParticleMechanicsApplication!\n"
            err_msg += "Available options are: \"residual_criterion\""
            raise Exception(err_msg)

        return convergence_criterion

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
                    KratosMultiphysics.Logger.PrintInfo('::[MPMSolver]:: ',\
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
        """Create the solution scheme for interested problem."""
        raise Exception("Solution Scheme creation must be implemented in the derived class.")

    def _create_solution_strategy(self):
        analysis_type = self.settings["analysis_type"].GetString()
        if analysis_type == "non_linear":
                solution_strategy = self._create_newton_raphson_strategy()
        else:
            err_msg =  "The requested analysis type \"" + analysis_type + "\" is not available!\n"
            err_msg += "Available options are: \"non_linear\""
            raise Exception(err_msg)
        return solution_strategy

    def _create_newton_raphson_strategy(self):
        computing_model_part = self.GetComputingModelPart()
        solution_scheme = self.get_solution_scheme()
        linear_solver = self.get_linear_solver()
        convergence_criterion = self.get_convergence_criterion()
        builder_and_solver = self.get_builder_and_solver()
        reform_dofs_at_each_step = False
        return KratosParticle.MPMResidualBasedNewtonRaphsonStrategy(computing_model_part,
                                                                        solution_scheme,
                                                                        linear_solver,
                                                                        convergence_criterion,
                                                                        builder_and_solver,
                                                                        self.settings["max_iteration"].GetInt(),
                                                                        self.settings["compute_reactions"].GetBool(),
                                                                        reform_dofs_at_each_step,
                                                                        self.settings["move_mesh_flag"].GetBool())

    def _set_buffer_size(self):
        current_buffer_size = self.grid_model_part.GetBufferSize()
        if self.min_buffer_size > current_buffer_size:
            self.grid_model_part.SetBufferSize(self.min_buffer_size)
        else:
            self.grid_model_part.SetBufferSize(current_buffer_size)

        current_buffer_size = self.initial_material_model_part.GetBufferSize()
        if self.min_buffer_size > current_buffer_size:
            self.initial_material_model_part.SetBufferSize(self.min_buffer_size)
        else:
            self.initial_material_model_part.SetBufferSize(current_buffer_size)

