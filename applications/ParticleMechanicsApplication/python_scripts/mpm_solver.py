from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics

# Import applications and dependencies
import KratosMultiphysics.ParticleMechanicsApplication as KratosParticle

# Importing the base class
from KratosMultiphysics.python_solver import PythonSolver

# Other imports
from KratosMultiphysics import auxiliary_solver_utilities

def CreateSolver(model, custom_settings):
    return MPMSolver(model, custom_settings)

class MPMSolver(PythonSolver):

    ### Solver constructor
    def __init__(self, model, custom_settings):
        self._validate_settings_in_baseclass=True # To be removed eventually
        super(MPMSolver, self).__init__(model, custom_settings)

        # Add model part containers
        self._AddModelPartContainers()

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
        self._AddVariablesToModelPart(self.grid_model_part)

        KratosMultiphysics.Logger.PrintInfo("::[MPMSolver]:: ", "Variables are added.")

    def ImportModelPart(self):
        # Read model part
        self._ModelPartReading()

        KratosMultiphysics.Logger.PrintInfo("::[MPMSolver]:: ","Models are imported.")

    def PrepareModelPart(self):
        # Set buffer size
        self._SetBufferSize()

        # Executes the check and prepare model process
        self.__ExecuteCheckAndPrepare()

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
        self._AddDofsToModelPart(self.grid_model_part)

        KratosMultiphysics.Logger.PrintInfo("::[MPMSolver]:: ","DOFs are added.")

    def Initialize(self):
        # The particle solution strategy is created here if it does not already exist.
        particle_solution_strategy = self._GetSolutionStrategy()
        particle_solution_strategy.SetEchoLevel(self.settings["echo_level"].GetInt())

        # Generate material points
        self._GenerateMaterialPoint()

        KratosMultiphysics.Logger.PrintInfo("::[MPMSolver]:: ","Solver is initialized correctly.")

    def AdvanceInTime(self, current_time):
        dt = self.__ComputeDeltaTime()
        new_time = current_time + dt
        self.grid_model_part.ProcessInfo[KratosMultiphysics.STEP] += 1
        self.grid_model_part.CloneTimeStep(new_time)

        return new_time

    def InitializeSolutionStep(self):
        self._SearchElement()
        self._GetSolutionStrategy().Initialize()
        self._GetSolutionStrategy().InitializeSolutionStep()

    def Predict(self):
        self._GetSolutionStrategy().Predict()

    def SolveSolutionStep(self):
        is_converged = self._GetSolutionStrategy().SolveSolutionStep()
        return is_converged

    def FinalizeSolutionStep(self):
        self._GetSolutionStrategy().FinalizeSolutionStep()

        self._GetSolutionStrategy().Clear()

    def Check(self):
        self._GetSolutionStrategy().Check()

    def Clear(self):
        self._GetSolutionStrategy().Clear()

    ### Solver special protected functions

    def _GetSolutionScheme(self):
        if not hasattr(self, '_solution_scheme'):
            self._solution_scheme = self._CreateSolutionScheme()
        return self._solution_scheme

    def _GetConvergenceCriteria(self):
        if not hasattr(self, '_convergence_criterion'):
            self._convergence_criterion = self._CreateConvergenceCriteria()
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
        if not hasattr(self, '_solution_strategy'):
            self._solution_strategy = self._CreateSolutionStrategy()
        return self._solution_strategy

    ### Solver protected functions

    def _GenerateMaterialPoint(self):
        pressure_dofs          = self.settings["pressure_dofs"].GetBool()
        axis_symmetric_flag    = self.settings["axis_symmetric_flag"].GetBool()

        # Assigning extra information to the main model part
        self.material_point_model_part.SetNodes(self.grid_model_part.GetNodes())
        self.material_point_model_part.ProcessInfo = self.grid_model_part.ProcessInfo
        self.material_point_model_part.SetBufferSize(self.grid_model_part.GetBufferSize())

        # Generate MP Element and Condition
        KratosParticle.GenerateMaterialPointElement(self.grid_model_part, self.initial_mesh_model_part, self.material_point_model_part, axis_symmetric_flag, pressure_dofs)
        KratosParticle.GenerateMaterialPointCondition(self.grid_model_part, self.initial_mesh_model_part, self.material_point_model_part)

    def _SearchElement(self):
        searching_alg_type = self.settings["element_search_settings"]["search_algorithm_type"].GetString()
        max_number_of_search_results = self.settings["element_search_settings"]["max_number_of_results"].GetInt()
        searching_tolerance          = self.settings["element_search_settings"]["searching_tolerance"].GetDouble()
        if (searching_alg_type == "bin_based"):
            KratosParticle.SearchElement(self.grid_model_part, self.material_point_model_part, max_number_of_search_results, searching_tolerance)
        else:
            err_msg  = "The requested searching algorithm \"" + searching_alg_type
            err_msg += "\" is not available for ParticleMechanicsApplication!\n"
            err_msg += "Available options are: \"bin_based\""
            raise Exception(err_msg)

    def _AddModelPartContainers(self):
        domain_size = self._GetDomainSize()
        if domain_size not in [2,3]:
            err_msg  = "The input \"domain_size\" is wrong!"
            err_msg += "Available options are: \"2\" or \"3\""
            raise Exception(err_msg)

        ## In MPM three model parts are needed
        # Material model part definition
        material_point_model_part_name = self.settings["model_part_name"].GetString()
        if not self.model.HasModelPart(material_point_model_part_name):
            self.material_point_model_part = self.model.CreateModelPart(material_point_model_part_name) # Equivalent to model_part3 in the old format
            self.material_point_model_part.ProcessInfo.SetValue(KratosMultiphysics.DOMAIN_SIZE, domain_size)

        # Initial material model part definition
        initial_mesh_model_part_name = "Initial_" + material_point_model_part_name
        if not self.model.HasModelPart(initial_mesh_model_part_name):
            self.initial_mesh_model_part = self.model.CreateModelPart(initial_mesh_model_part_name) #Equivalent to model_part2 in the old format
            self.initial_mesh_model_part.ProcessInfo.SetValue(KratosMultiphysics.DOMAIN_SIZE, domain_size)

        # Grid model part definition
        if not self.model.HasModelPart("Background_Grid"):
            self.grid_model_part = self.model.CreateModelPart("Background_Grid") #Equivalent to model_part1 in the old format
            self.grid_model_part.ProcessInfo.SetValue(KratosMultiphysics.DOMAIN_SIZE, domain_size)

    def _AddVariablesToModelPart(self, model_part):
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
        auxiliary_solver_utilities.AddVariables(model_part, self.settings["auxiliary_variables_list"])

        # Add variables for specific cases
        if self.settings["pressure_dofs"].GetBool():
            # add specific variables for the problem (pressure dofs)
            model_part.AddNodalSolutionStepVariable(KratosParticle.PRESSURE_REACTION)
            model_part.AddNodalSolutionStepVariable(KratosParticle.NODAL_MPRESSURE)

    def _AddDynamicVariables(self, model_part):
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VELOCITY)
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.ACCELERATION)

    def _ModelPartReading(self):
        # reading the model part of the background grid
        if(self.settings["grid_model_import_settings"]["input_type"].GetString() == "mdpa"):
            self._ImportModelPart(self.grid_model_part, self.settings["grid_model_import_settings"])
        else:
            raise Exception("Other input options are not implemented yet.")

        # reading the model part of the material point
        if(self.settings["model_import_settings"]["input_type"].GetString() == "mdpa"):
            self._ImportModelPart(self.initial_mesh_model_part, self.settings["model_import_settings"])
        else:
            raise Exception("Other input options are not implemented yet.")

    def _AddDofsToModelPart(self, model_part):
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_X, KratosMultiphysics.REACTION_X, model_part)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_Y, KratosMultiphysics.REACTION_Y, model_part)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_Z, KratosMultiphysics.REACTION_Z, model_part)

        if self.settings["pressure_dofs"].GetBool():
            KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.PRESSURE, KratosParticle.PRESSURE_REACTION, model_part)

        # Add dofs that the user defined in the ProjectParameters
        auxiliary_solver_utilities.AddDofs(model_part, self.settings["auxiliary_dofs_list"], self.settings["auxiliary_reaction_list"])

    def _GetDomainSize(self):
        if not hasattr(self, '_domain_size'):
            self._domain_size = self.settings["domain_size"].GetInt()
        return self._domain_size

    def _GetConvergenceCriteriaSettings(self):
        # Create an auxiliary Kratos parameters object to store the convergence settings.
        conv_params = KratosMultiphysics.Parameters("{}")
        conv_params.AddValue("convergence_criterion",self.settings["convergence_criterion"])
        conv_params.AddValue("echo_level",self.settings["echo_level"])
        conv_params.AddValue("displacement_relative_tolerance",self.settings["displacement_relative_tolerance"])
        conv_params.AddValue("displacement_absolute_tolerance",self.settings["displacement_absolute_tolerance"])
        conv_params.AddValue("residual_relative_tolerance",self.settings["residual_relative_tolerance"])
        conv_params.AddValue("residual_absolute_tolerance",self.settings["residual_absolute_tolerance"])

        return conv_params

    def _CreateConvergenceCriteria(self):
        convergence_criterion_parameters = self._GetConvergenceCriteriaSettings()
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

    def _CreateLinearSolver(self):
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

    def _CreateBuilderAndSolver(self):
        linear_solver = self._GetLinearSolver()
        if self.settings["block_builder"].GetBool():
            builder_and_solver = KratosMultiphysics.ResidualBasedBlockBuilderAndSolver(linear_solver)
        else:
            builder_and_solver = KratosMultiphysics.ResidualBasedEliminationBuilderAndSolver(linear_solver)
        return builder_and_solver

    def _CreateSolutionScheme(self):
        """Create the solution scheme for interested problem."""
        raise Exception("Solution Scheme creation must be implemented in the derived class.")

    def _CreateSolutionStrategy(self):
        analysis_type = self.settings["analysis_type"].GetString()
        if analysis_type == "non_linear":
                solution_strategy = self._CreateNewtonRaphsonStrategy()
        elif analysis_type == 'linear':
                self.material_point_model_part.ProcessInfo.SetValue(KratosParticle.IGNORE_GEOMETRIC_STIFFNESS, True)
                solution_strategy = self._CreateLinearStrategy();
        else:
            err_msg =  "The requested analysis type \"" + analysis_type + "\" is not available!\n"
            err_msg += "Available options are: \"linear\", \"non_linear\""
            raise Exception(err_msg)
        return solution_strategy

    def _CreateNewtonRaphsonStrategy(self):
        computing_model_part = self.GetComputingModelPart()
        solution_scheme = self._GetSolutionScheme()
        linear_solver = self._GetLinearSolver()
        convergence_criterion = self._GetConvergenceCriteria()
        builder_and_solver = self._GetBuilderAndSolver()
        reform_dofs_at_each_step = False ## hard-coded, but can be changed upon implementation
        return KratosParticle.MPMResidualBasedNewtonRaphsonStrategy(computing_model_part,
                                                                        solution_scheme,
                                                                        linear_solver,
                                                                        convergence_criterion,
                                                                        builder_and_solver,
                                                                        self.settings["max_iteration"].GetInt(),
                                                                        self.settings["compute_reactions"].GetBool(),
                                                                        reform_dofs_at_each_step,
                                                                        self.settings["move_mesh_flag"].GetBool())

    def _CreateLinearStrategy(self):
        computing_model_part = self.GetComputingModelPart()
        solution_scheme = self._GetSolutionScheme()
        linear_solver = self._GetLinearSolver()
        reform_dofs_at_each_step = False ## hard-coded, but can be changed upon implementation
        calc_norm_dx_flag = False ## hard-coded, but can be changed upon implementation
        return KratosMultiphysics.ResidualBasedLinearStrategy(computing_model_part,
                                                              solution_scheme,
                                                              linear_solver,
                                                              self.settings["compute_reactions"].GetBool(),
                                                              reform_dofs_at_each_step,
                                                              calc_norm_dx_flag,
                                                              self.settings["move_mesh_flag"].GetBool())

    def _SetBufferSize(self):
        current_buffer_size = self.grid_model_part.GetBufferSize()
        if self.min_buffer_size > current_buffer_size:
            self.grid_model_part.SetBufferSize(self.min_buffer_size)
        else:
            self.grid_model_part.SetBufferSize(current_buffer_size)

        current_buffer_size = self.initial_mesh_model_part.GetBufferSize()
        if self.min_buffer_size > current_buffer_size:
            self.initial_mesh_model_part.SetBufferSize(self.min_buffer_size)
        else:
            self.initial_mesh_model_part.SetBufferSize(current_buffer_size)

    ### Solver private functions

    def __ComputeDeltaTime(self):
        return self.settings["time_stepping"]["time_step"].GetDouble()

    def __ExecuteCheckAndPrepare(self):
        # Specific active node and element check for particle MPM solver
        for node in self.grid_model_part.Nodes:
            if (node.Is(KratosMultiphysics.ACTIVE)):
                KratosMultiphysics.Logger.PrintInfo("::[MPMSolver]:: ","WARNING: This grid node have been set active: ", node.Id)

        # Setting active initial elements
        KratosMultiphysics.VariableUtils().SetFlag(KratosMultiphysics.ACTIVE, True, self.initial_mesh_model_part.Elements)

        # Read material property
        materials_imported = self.__ImportConstitutiveLaws()
        if materials_imported:
            KratosMultiphysics.Logger.PrintInfo("::[MPMSolver]:: ","Constitutive law was successfully imported.")
        else:
            KratosMultiphysics.Logger.PrintWarning("::[MPMSolver]:: ","Constitutive law was not imported.")

        # Clone property of model_part2 to model_part3
        self.material_point_model_part.Properties = self.initial_mesh_model_part.Properties

    def __ImportConstitutiveLaws(self):
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

