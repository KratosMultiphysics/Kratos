from __future__ import print_function, absolute_import, division  # makes KM backward compatible with python 2.6 and 2.7
#import kratos core and applications
import KratosMultiphysics as KM

# Check that applications were imported in the main script
KM.CheckRegisteredApplications("StructuralMechanicsApplication")
KM.CheckRegisteredApplications("ContactStructuralMechanicsApplication")

# Import applications
import KratosMultiphysics.StructuralMechanicsApplication as SMA
import KratosMultiphysics.ContactStructuralMechanicsApplication as CSMA

# Import the implicit solver (the explicit one is derived from it)
import structural_mechanics_implicit_dynamic_solver

def CreateSolver(model, custom_settings):
    return ContactImplicitMechanicalSolver(model, custom_settings)

class ContactImplicitMechanicalSolver(structural_mechanics_implicit_dynamic_solver.ImplicitMechanicalSolver):
    """The structural mechanics contact implicit dynamic solver.

    This class creates the mechanical solvers for contact implicit dynamic analysis.
    It currently supports Newmark, Bossak and dynamic relaxation schemes.

    Public member variables:
    dynamic_settings -- settings for the implicit dynamic solvers.

    See structural_mechanics_solver.py for more information.
    """
    def __init__(self, model, custom_settings):

        ##settings string in json format
        contact_settings = KM.Parameters("""
        {
            "contact_settings" :
            {
                "mortar_type"                                       : "",
                "condn_convergence_criterion"                       : false,
                "fancy_convergence_criterion"                       : true,
                "print_convergence_criterion"                       : false,
                "ensure_contact"                                    : false,
                "frictional_decomposed"                             : true,
                "gidio_debug"                                       : false,
                "adaptative_strategy"                               : false,
                "split_factor"                                      : 10.0,
                "max_number_splits"                                 : 3,
                "inner_loop_iterations"                             : 5,
                "contact_displacement_relative_tolerance"           : 1.0e-4,
                "contact_displacement_absolute_tolerance"           : 1.0e-9,
                "contact_residual_relative_tolerance"               : 1.0e-4,
                "contact_residual_absolute_tolerance"               : 1.0e-9,
                "frictional_contact_displacement_relative_tolerance": 1.0e-4,
                "frictional_contact_displacement_absolute_tolerance": 1.0e-9,
                "frictional_contact_residual_relative_tolerance"    : 1.0e-4,
                "frictional_contact_residual_absolute_tolerance"    : 1.0e-9,
                "simplified_semi_smooth_newton"                     : false,
                "rescale_linear_solver"                             : false,
                "use_mixed_ulm_solver"                              : true,
                "mixed_ulm_solver_parameters" :
                {
                    "solver_type"          : "mixed_ulm_linear_solver",
                    "tolerance"            : 1.0e-6,
                    "max_iteration_number" : 200,
                    "echo_level"           : 0
                }
            }
        }
        """)

        ## Overwrite the default settings with user-provided parameters
        self.settings = custom_settings
        self.validate_and_transfer_matching_settings(self.settings, contact_settings)
        self.contact_settings = contact_settings["contact_settings"]

        # Linear solver settings
        if (self.settings.Has("linear_solver_settings")):
            self.linear_solver_settings = self.settings["linear_solver_settings"]
        else:
            self.linear_solver_settings = KM.Parameters("""{}""")

        # Construct the base solver.
        super(ContactImplicitMechanicalSolver, self).__init__(model, self.settings)

        # Setting default configurations true by default
        if (self.settings["clear_storage"].GetBool() == False):
            self.print_on_rank_zero("Clear storage", "Storage must be cleared each step. Switching to True")
            self.settings["clear_storage"].SetBool(True)
        if (self.settings["reform_dofs_at_each_step"].GetBool() == False):
            self.print_on_rank_zero("Reform DoFs", "DoF must be reformed each time step. Switching to True")
            self.settings["reform_dofs_at_each_step"].SetBool(True)
        if (self.settings["block_builder"].GetBool() == False):
            self.print_on_rank_zero("Builder and solver", "EliminationBuilderAndSolver can not used with the current implementation. Switching to BlockBuilderAndSolver")
            self.settings["block_builder"].SetBool(True)

        # Setting echo level
        self.echo_level =  self.settings["echo_level"].GetInt()

        # Initialize the processes list
        self.processes_list = None

        # Initialize the post process
        self.post_process = None

        self.print_on_rank_zero("::[Contact Mechanical Implicit Dynamic Solver]:: ", "Construction of ContactMechanicalSolver finished")

    def AddVariables(self):

        super(ContactImplicitMechanicalSolver, self).AddVariables()

        mortar_type = self.contact_settings["mortar_type"].GetString()
        if mortar_type != "":
            self.main_model_part.AddNodalSolutionStepVariable(KM.NORMAL)  # Add normal
            self.main_model_part.AddNodalSolutionStepVariable(KM.NODAL_H) # Add nodal size variable
            if mortar_type == "ALMContactFrictionless":
                self.main_model_part.AddNodalSolutionStepVariable(CSMA.LAGRANGE_MULTIPLIER_CONTACT_PRESSURE) # Add normal contact stress
                self.main_model_part.AddNodalSolutionStepVariable(CSMA.WEIGHTED_GAP)                         # Add normal contact gap
                self.main_model_part.AddNodalSolutionStepVariable(CSMA.WEIGHTED_SCALAR_RESIDUAL)             # Add scalar LM residual
            elif mortar_type == "ALMContactFrictionlessComponents":
                self.main_model_part.AddNodalSolutionStepVariable(KM.VECTOR_LAGRANGE_MULTIPLIER)             # Add normal contact stress
                self.main_model_part.AddNodalSolutionStepVariable(CSMA.WEIGHTED_GAP)                         # Add normal contact gap
                self.main_model_part.AddNodalSolutionStepVariable(CSMA.WEIGHTED_VECTOR_RESIDUAL)             # Add vector LM residual
            elif mortar_type == "ALMContactFrictional":
                self.main_model_part.AddNodalSolutionStepVariable(KM.VECTOR_LAGRANGE_MULTIPLIER)             # Add normal contact stress
                self.main_model_part.AddNodalSolutionStepVariable(CSMA.WEIGHTED_GAP)                         # Add normal contact gap
                self.main_model_part.AddNodalSolutionStepVariable(CSMA.WEIGHTED_SLIP)                        # Add contact slip
                self.main_model_part.AddNodalSolutionStepVariable(CSMA.WEIGHTED_VECTOR_RESIDUAL)             # Add vector LM residual
            elif mortar_type == "ScalarMeshTying":
                self.main_model_part.AddNodalSolutionStepVariable(KM.SCALAR_LAGRANGE_MULTIPLIER)             # Add scalar LM
                self.main_model_part.AddNodalSolutionStepVariable(CSMA.WEIGHTED_SCALAR_RESIDUAL)             # Add scalar LM residual
            elif mortar_type == "ComponentsMeshTying":
                self.main_model_part.AddNodalSolutionStepVariable(KM.VECTOR_LAGRANGE_MULTIPLIER)             # Add vector LM
                self.main_model_part.AddNodalSolutionStepVariable(CSMA.WEIGHTED_VECTOR_RESIDUAL)             # Add vector LM residual

        self.print_on_rank_zero("::[Contact Mechanical Implicit Dynamic Solver]:: ", "Variables ADDED")

    def AddDofs(self):

        super(ContactImplicitMechanicalSolver, self).AddDofs()

        mortar_type = self.contact_settings["mortar_type"].GetString()
        if (mortar_type == "ALMContactFrictionless"):                                                      # TODO Remove WEIGHTED_SCALAR_RESIDUAL in case of check for reaction is defined
            KM.VariableUtils().AddDof(CSMA.LAGRANGE_MULTIPLIER_CONTACT_PRESSURE, CSMA.WEIGHTED_SCALAR_RESIDUAL, self.main_model_part)
        elif (mortar_type == "ALMContactFrictional" or mortar_type == "ALMContactFrictionlessComponents"): # TODO Remove WEIGHTED_VECTOR_RESIDUAL in case of check for reaction is defined
            KM.VariableUtils().AddDof(KM.VECTOR_LAGRANGE_MULTIPLIER_X, CSMA.WEIGHTED_VECTOR_RESIDUAL_X, self.main_model_part)
            KM.VariableUtils().AddDof(KM.VECTOR_LAGRANGE_MULTIPLIER_Y, CSMA.WEIGHTED_VECTOR_RESIDUAL_Y, self.main_model_part)
            KM.VariableUtils().AddDof(KM.VECTOR_LAGRANGE_MULTIPLIER_Z, CSMA.WEIGHTED_VECTOR_RESIDUAL_Z, self.main_model_part)
        elif (mortar_type == "ScalarMeshTying"):
            KM.VariableUtils().AddDof(KM.SCALAR_LAGRANGE_MULTIPLIER,CSMA.WEIGHTED_SCALAR_RESIDUAL, self.main_model_part)
        elif (mortar_type == "ComponentsMeshTying"):
            KM.VariableUtils().AddDof(KM.VECTOR_LAGRANGE_MULTIPLIER_X, CSMA.WEIGHTED_VECTOR_RESIDUAL_X, self.main_model_part)
            KM.VariableUtils().AddDof(KM.VECTOR_LAGRANGE_MULTIPLIER_Y, CSMA.WEIGHTED_VECTOR_RESIDUAL_Y, self.main_model_part)
            KM.VariableUtils().AddDof(KM.VECTOR_LAGRANGE_MULTIPLIER_Z, CSMA.WEIGHTED_VECTOR_RESIDUAL_Z, self.main_model_part)

        self.print_on_rank_zero("::[Contact Mechanical Implicit Dynamic Solver]:: ", "DOF's ADDED")

    def Initialize(self):
        super(ContactImplicitMechanicalSolver, self).Initialize() # The mechanical solver is created here.

        # We set the flag INTERACTION
        computing_model_part = self.GetComputingModelPart()
        if self.contact_settings["simplified_semi_smooth_newton"].GetBool() is True:
            computing_model_part.Set(KM.INTERACTION, False)
        else:
            computing_model_part.Set(KM.INTERACTION, True)

    def Solve(self):
        if self.settings["clear_storage"].GetBool():
            self.Clear()

        mechanical_solution_strategy = self.get_mechanical_solution_strategy()

        # The steps of the solve are Initialize(), InitializeSolutionStep(), Predict(), SolveSolutionStep(), FinalizeSolutionStep()
        mechanical_solution_strategy.Solve()
        #mechanical_solution_strategy.Initialize()
        #mechanical_solution_strategy.InitializeSolutionStep()
        #mechanical_solution_strategy.Predict()
        #mechanical_solution_strategy.SolveSolutionStep()
        #mechanical_solution_strategy.FinalizeSolutionStep()

    def AddProcessesList(self, processes_list):
        self.processes_list = CSMA.ProcessFactoryUtility(processes_list)

    def AddPostProcess(self, post_process):
        self.post_process = CSMA.ProcessFactoryUtility(post_process)

    def print_on_rank_zero(self, *args):
        # This function will be overridden in the trilinos-solvers
        KM.Logger.PrintInfo(" ".join(map(str,args)))

    def print_warning_on_rank_zero(self, *args):
        # This function will be overridden in the trilinos-solvers
        KM.Logger.PrintWarning(" ".join(map(str,args)))

    def _create_linear_solver(self):
        linear_solver = super(ContactImplicitMechanicalSolver, self)._create_linear_solver()
        if (self.contact_settings["rescale_linear_solver"].GetBool() is True):
            linear_solver = KM.ScalingSolver(linear_solver, False)
        mortar_type = self.contact_settings["mortar_type"].GetString()
        if (mortar_type == "ALMContactFrictional" or mortar_type == "ALMContactFrictionlessComponents"):
            if (self.contact_settings["use_mixed_ulm_solver"].GetBool() == True):
                self.print_on_rank_zero("::[Contact Mechanical Implicit Dynamic Solver]:: ", "Using MixedULMLinearSolver, definition of ALM parameters recommended")
                name_mixed_solver = self.contact_settings["mixed_ulm_solver_parameters"]["solver_type"].GetString()
                if (name_mixed_solver == "mixed_ulm_linear_solver"):
                    linear_solver_name = self.settings["linear_solver_settings"]["solver_type"].GetString()
                    if (linear_solver_name == "AMGCL" or linear_solver_name == "AMGCLSolver"):
                        amgcl_param = KM.Parameters("""
                        {
                            "solver_type"                    : "AMGCL",
                            "smoother_type"                  : "ilu0",
                            "krylov_type"                    : "lgmres",
                            "coarsening_type"                : "aggregation",
                            "max_iteration"                  : 100,
                            "provide_coordinates"            : false,
                            "gmres_krylov_space_dimension"   : 100,
                            "verbosity"                      : 1,
                            "tolerance"                      : 1e-6,
                            "scaling"                        : false,
                            "block_size"                     : 3,
                            "use_block_matrices_if_possible" : true,
                            "coarse_enough"                  : 500
                        }
                        """)
                        amgcl_param["block_size"].SetInt(self.main_model_part.ProcessInfo[KM.DOMAIN_SIZE])
                        self.linear_solver_settings.RecursivelyValidateAndAssignDefaults(amgcl_param)
                        linear_solver = KM.AMGCLSolver(self.linear_solver_settings)
                    mixed_ulm_solver = CSMA.MixedULMLinearSolver(linear_solver, self.contact_settings["mixed_ulm_solver_parameters"])
                    return mixed_ulm_solver
                else:
                    self.print_on_rank_zero("::[Contact Mechanical Implicit Dynamic Solver]:: ", "Mixed solver not available: " + name_mixed_solver + ". Using not mixed linear solver")
                    return linear_solver
            else:
                return linear_solver
        else:
            return linear_solver

    def _get_convergence_criterion_settings(self):
        # Create an auxiliary Kratos parameters object to store the convergence settings.
        if (self.contact_settings["fancy_convergence_criterion"].GetBool() is True):
            table = KM.TableStreamUtility()
            table.SetOnProcessInfo(self.main_model_part.ProcessInfo)

        conv_params = KM.Parameters("{}")
        conv_params.AddValue("convergence_criterion", self.settings["convergence_criterion"])
        conv_params.AddValue("rotation_dofs", self.settings["rotation_dofs"])
        conv_params.AddValue("echo_level", self.settings["echo_level"])
        conv_params.AddValue("displacement_relative_tolerance", self.settings["displacement_relative_tolerance"])
        conv_params.AddValue("displacement_absolute_tolerance", self.settings["displacement_absolute_tolerance"])
        conv_params.AddValue("residual_relative_tolerance", self.settings["residual_relative_tolerance"])
        conv_params.AddValue("residual_absolute_tolerance", self.settings["residual_absolute_tolerance"])
        conv_params.AddValue("contact_displacement_relative_tolerance", self.contact_settings["contact_displacement_relative_tolerance"])
        conv_params.AddValue("contact_displacement_absolute_tolerance", self.contact_settings["contact_displacement_absolute_tolerance"])
        conv_params.AddValue("contact_residual_relative_tolerance", self.contact_settings["contact_residual_relative_tolerance"])
        conv_params.AddValue("contact_residual_absolute_tolerance", self.contact_settings["contact_residual_absolute_tolerance"])
        conv_params.AddValue("frictional_contact_displacement_relative_tolerance", self.contact_settings["frictional_contact_displacement_relative_tolerance"])
        conv_params.AddValue("frictional_contact_displacement_absolute_tolerance", self.contact_settings["frictional_contact_displacement_absolute_tolerance"])
        conv_params.AddValue("frictional_contact_residual_relative_tolerance", self.contact_settings["frictional_contact_residual_relative_tolerance"])
        conv_params.AddValue("frictional_contact_residual_absolute_tolerance", self.contact_settings["frictional_contact_residual_absolute_tolerance"])
        conv_params.AddValue("mortar_type", self.contact_settings["mortar_type"])
        conv_params.AddValue("condn_convergence_criterion", self.contact_settings["condn_convergence_criterion"])
        conv_params.AddValue("print_convergence_criterion", self.contact_settings["print_convergence_criterion"])
        conv_params.AddValue("ensure_contact", self.contact_settings["ensure_contact"])
        conv_params.AddValue("frictional_decomposed", self.contact_settings["frictional_decomposed"])
        conv_params.AddValue("gidio_debug", self.contact_settings["gidio_debug"])

        return conv_params

    def _create_convergence_criterion(self):
        import contact_convergence_criteria_factory
        convergence_criterion = contact_convergence_criteria_factory.convergence_criterion(self._get_convergence_criterion_settings())
        return convergence_criterion.mechanical_convergence_criterion

    def _create_builder_and_solver(self):
        if self.contact_settings["mortar_type"].GetString() != "":
            linear_solver = self.get_linear_solver()
            if self.settings["block_builder"].GetBool():
                if self.settings["multi_point_constraints_used"].GetBool():
                    builder_and_solver = CSMA.ContactResidualBasedBlockBuilderAndSolverWithConstraints(linear_solver)
                else:
                    builder_and_solver = CSMA.ContactResidualBasedBlockBuilderAndSolver(linear_solver)
            else:
                raise Exception("Contact not compatible with EliminationBuilderAndSolver")
        else:
            builder_and_solver = super(ContactImplicitMechanicalSolver, self)._create_builder_and_solver()

        return builder_and_solver

    def _create_mechanical_solution_strategy(self):
        if self.contact_settings["mortar_type"].GetString() != "":
            if self.settings["analysis_type"].GetString() == "linear":
                mechanical_solution_strategy = self._create_linear_strategy()
            else:
                if(self.settings["line_search"].GetBool()):
                    mechanical_solution_strategy = self._create_contact_line_search_strategy()
                else:
                    mechanical_solution_strategy = self._create_contact_newton_raphson_strategy()
        else:
            mechanical_solution_strategy = super(ContactImplicitMechanicalSolver, self)._create_mechanical_solution_strategy()

        return mechanical_solution_strategy

    def _create_contact_line_search_strategy(self):
        computing_model_part = self.GetComputingModelPart()
        self.mechanical_scheme = self.get_solution_scheme()
        self.linear_solver = self.get_linear_solver()
        self.mechanical_convergence_criterion = self.get_convergence_criterion()
        self.builder_and_solver = self.get_builder_and_solver()
        newton_parameters = KM.Parameters("""{}""")
        return CSMA.LineSearchContactStrategy(computing_model_part,
                                                self.mechanical_scheme,
                                                self.linear_solver,
                                                self.mechanical_convergence_criterion,
                                                self.builder_and_solver,
                                                self.settings["max_iteration"].GetInt(),
                                                self.settings["compute_reactions"].GetBool(),
                                                self.settings["reform_dofs_at_each_step"].GetBool(),
                                                self.settings["move_mesh_flag"].GetBool(),
                                                newton_parameters
                                                )
    def _create_contact_newton_raphson_strategy(self):
        computing_model_part = self.GetComputingModelPart()
        self.mechanical_scheme = self.get_solution_scheme()
        self.linear_solver = self.get_linear_solver()
        self.mechanical_convergence_criterion = self.get_convergence_criterion()
        self.builder_and_solver = self.get_builder_and_solver()
        newton_parameters = KM.Parameters("""{}""")
        newton_parameters.AddValue("adaptative_strategy", self.contact_settings["adaptative_strategy"])
        newton_parameters.AddValue("split_factor", self.contact_settings["split_factor"])
        newton_parameters.AddValue("max_number_splits", self.contact_settings["max_number_splits"])
        newton_parameters.AddValue("inner_loop_iterations", self.contact_settings["inner_loop_iterations"])
        return CSMA.ResidualBasedNewtonRaphsonContactStrategy(computing_model_part,
                                                                self.mechanical_scheme,
                                                                self.linear_solver,
                                                                self.mechanical_convergence_criterion,
                                                                self.builder_and_solver,
                                                                self.settings["max_iteration"].GetInt(),
                                                                self.settings["compute_reactions"].GetBool(),
                                                                self.settings["reform_dofs_at_each_step"].GetBool(),
                                                                self.settings["move_mesh_flag"].GetBool(),
                                                                newton_parameters,
                                                                self.processes_list,
                                                                self.post_process
                                                                )
