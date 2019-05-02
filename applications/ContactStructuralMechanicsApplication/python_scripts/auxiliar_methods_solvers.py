from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# Importing the Kratos Library
import KratosMultiphysics as KM

import KratosMultiphysics.StructuralMechanicsApplication as SMA
import KratosMultiphysics.ContactStructuralMechanicsApplication as CSMA


def print_on_rank_zero(*args):
    # This function will be overridden in the trilinos-solvers
    KM.Logger.PrintInfo(" ".join(map(str, args)))


def print_warning_on_rank_zero(*args):
    # This function will be overridden in the trilinos-solvers
    KM.Logger.PrintWarning(" ".join(map(str, args)))


def AuxiliarContactSettings():
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
            "compute_dynamic_factor"                            : false,
            "gidio_debug"                                       : false,
            "adaptative_strategy"                               : false,
            "split_factor"                                      : 10.0,
            "max_number_splits"                                 : 3,
            "inner_loop_iterations"                             : 5,
            "inner_loop_adaptive"                               : false,
            "contact_displacement_relative_tolerance"           : 1.0e-4,
            "contact_displacement_absolute_tolerance"           : 1.0e-9,
            "contact_residual_relative_tolerance"               : 1.0e-4,
            "contact_residual_absolute_tolerance"               : 1.0e-9,
            "frictional_contact_displacement_relative_tolerance": 1.0e-4,
            "frictional_contact_displacement_absolute_tolerance": 1.0e-9,
            "frictional_contact_residual_relative_tolerance"    : 1.0e-4,
            "frictional_contact_residual_absolute_tolerance"    : 1.0e-9,
            "silent_strategy"                                   : true,
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

    return contact_settings


def AuxiliarExplicitContactSettings():
    contact_settings = KM.Parameters("""
    {
        "contact_settings" :
        {
            "mortar_type"                                       : "",
            "compute_dynamic_factor"                            : true,
            "ensure_contact"                                    : false,
            "silent_strategy"                                   : false,
            "delta_time_factor_for_contact"                     : 5.0e-1
        }
    }
    """)

    return contact_settings


def AuxiliarSetSettings(settings, contact_settings):
    if not settings["clear_storage"].GetBool():
        print_on_rank_zero("Clear storage", "Storage must be cleared each step. Switching to True")
        settings["clear_storage"].SetBool(True)
    if not settings["reform_dofs_at_each_step"].GetBool():
        print_on_rank_zero("Reform DoFs", "DoF must be reformed each time step. Switching to True")
        settings["reform_dofs_at_each_step"].SetBool(True)
    if not settings["block_builder"].GetBool():
        print_on_rank_zero("Builder and solver", "EliminationBuilderAndSolver can not used with the current implementation. Switching to BlockBuilderAndSolver")
        settings["block_builder"].SetBool(True)

    return settings


def AuxiliarAddVariables(main_model_part, mortar_type = ""):
    if mortar_type != "":
        main_model_part.AddNodalSolutionStepVariable(KM.NORMAL)  # Add normal
        main_model_part.AddNodalSolutionStepVariable(KM.NODAL_H) # Add nodal size variable
        if mortar_type == "PenaltyContactFrictionless":
            main_model_part.AddNodalSolutionStepVariable(CSMA.WEIGHTED_GAP)                         # Add normal contact gap
        elif mortar_type == "PenaltyContactFrictional":
            main_model_part.AddNodalSolutionStepVariable(CSMA.WEIGHTED_GAP)                         # Add normal contact gap
            main_model_part.AddNodalSolutionStepVariable(CSMA.WEIGHTED_SLIP)                        # Add contact slip
        elif mortar_type == "ALMContactFrictionless":
            main_model_part.AddNodalSolutionStepVariable(CSMA.LAGRANGE_MULTIPLIER_CONTACT_PRESSURE) # Add normal contact stress
            main_model_part.AddNodalSolutionStepVariable(CSMA.WEIGHTED_GAP)                         # Add normal contact gap
            main_model_part.AddNodalSolutionStepVariable(CSMA.WEIGHTED_SCALAR_RESIDUAL)             # Add scalar LM residual
        elif mortar_type == "ALMContactFrictionlessComponents":
            main_model_part.AddNodalSolutionStepVariable(KM.VECTOR_LAGRANGE_MULTIPLIER)             # Add normal contact stress
            main_model_part.AddNodalSolutionStepVariable(CSMA.WEIGHTED_GAP)                         # Add normal contact gap
            main_model_part.AddNodalSolutionStepVariable(CSMA.WEIGHTED_VECTOR_RESIDUAL)             # Add vector LM residual
        elif mortar_type == "ALMContactFrictional":
            main_model_part.AddNodalSolutionStepVariable(KM.VECTOR_LAGRANGE_MULTIPLIER)             # Add normal contact stress
            main_model_part.AddNodalSolutionStepVariable(CSMA.WEIGHTED_GAP)                         # Add normal contact gap
            main_model_part.AddNodalSolutionStepVariable(CSMA.WEIGHTED_SLIP)                        # Add contact slip
            main_model_part.AddNodalSolutionStepVariable(CSMA.WEIGHTED_VECTOR_RESIDUAL)             # Add vector LM residual
        elif mortar_type == "ScalarMeshTying":
            main_model_part.AddNodalSolutionStepVariable(KM.SCALAR_LAGRANGE_MULTIPLIER)             # Add scalar LM
            main_model_part.AddNodalSolutionStepVariable(CSMA.WEIGHTED_SCALAR_RESIDUAL)             # Add scalar LM residual
        elif mortar_type == "ComponentsMeshTying":
            main_model_part.AddNodalSolutionStepVariable(KM.VECTOR_LAGRANGE_MULTIPLIER)             # Add vector LM
            main_model_part.AddNodalSolutionStepVariable(CSMA.WEIGHTED_VECTOR_RESIDUAL)             # Add vector LM residual


def AuxiliarAddDofs(main_model_part, mortar_type = ""):
    if mortar_type == "ALMContactFrictionless":                                                      # TODO Remove WEIGHTED_SCALAR_RESIDUAL in case of check for reaction is defined
        KM.VariableUtils().AddDof(CSMA.LAGRANGE_MULTIPLIER_CONTACT_PRESSURE, CSMA.WEIGHTED_SCALAR_RESIDUAL, main_model_part)
    elif mortar_type == "ALMContactFrictional" or mortar_type == "ALMContactFrictionlessComponents": # TODO Remove WEIGHTED_VECTOR_RESIDUAL in case of check for reaction is defined
        KM.VariableUtils().AddDof(KM.VECTOR_LAGRANGE_MULTIPLIER_X, CSMA.WEIGHTED_VECTOR_RESIDUAL_X, main_model_part)
        KM.VariableUtils().AddDof(KM.VECTOR_LAGRANGE_MULTIPLIER_Y, CSMA.WEIGHTED_VECTOR_RESIDUAL_Y, main_model_part)
        KM.VariableUtils().AddDof(KM.VECTOR_LAGRANGE_MULTIPLIER_Z, CSMA.WEIGHTED_VECTOR_RESIDUAL_Z, main_model_part)
    elif mortar_type == "ScalarMeshTying":
        KM.VariableUtils().AddDof(KM.SCALAR_LAGRANGE_MULTIPLIER,CSMA.WEIGHTED_SCALAR_RESIDUAL, main_model_part)
    elif mortar_type == "ComponentsMeshTying":
        KM.VariableUtils().AddDof(KM.VECTOR_LAGRANGE_MULTIPLIER_X, CSMA.WEIGHTED_VECTOR_RESIDUAL_X, main_model_part)
        KM.VariableUtils().AddDof(KM.VECTOR_LAGRANGE_MULTIPLIER_Y, CSMA.WEIGHTED_VECTOR_RESIDUAL_Y, main_model_part)
        KM.VariableUtils().AddDof(KM.VECTOR_LAGRANGE_MULTIPLIER_Z, CSMA.WEIGHTED_VECTOR_RESIDUAL_Z, main_model_part)


def AuxiliarSolve(mechanical_solution_strategy):
    # The steps of the solve are Initialize(), InitializeSolutionStep(), Predict(), SolveSolutionStep(), FinalizeSolutionStep()
    mechanical_solution_strategy.Solve()
    # mechanical_solution_strategy.Initialize()
    # mechanical_solution_strategy.InitializeSolutionStep()
    # mechanical_solution_strategy.Predict()
    # mechanical_solution_strategy.SolveSolutionStep()
    # mechanical_solution_strategy.FinalizeSolutionStep()


def AuxiliarCreateConvergenceParameters(main_model_part, settings, contact_settings):
    # Create an auxiliary Kratos parameters object to store the convergence settings.
        if contact_settings["fancy_convergence_criterion"].GetBool():
            table = KM.TableStreamUtility()
            table.SetOnProcessInfo(main_model_part.ProcessInfo)

        conv_params = KM.Parameters("{}")
        conv_params.AddValue("convergence_criterion", settings["convergence_criterion"])
        conv_params.AddValue("rotation_dofs", settings["rotation_dofs"])
        conv_params.AddValue("echo_level", settings["echo_level"])
        conv_params.AddValue("displacement_relative_tolerance", settings["displacement_relative_tolerance"])
        conv_params.AddValue("displacement_absolute_tolerance", settings["displacement_absolute_tolerance"])
        conv_params.AddValue("residual_relative_tolerance", settings["residual_relative_tolerance"])
        conv_params.AddValue("residual_absolute_tolerance", settings["residual_absolute_tolerance"])
        conv_params.AddValue("contact_displacement_relative_tolerance", contact_settings["contact_displacement_relative_tolerance"])
        conv_params.AddValue("contact_displacement_absolute_tolerance", contact_settings["contact_displacement_absolute_tolerance"])
        conv_params.AddValue("contact_residual_relative_tolerance", contact_settings["contact_residual_relative_tolerance"])
        conv_params.AddValue("contact_residual_absolute_tolerance", contact_settings["contact_residual_absolute_tolerance"])
        conv_params.AddValue("frictional_contact_displacement_relative_tolerance", contact_settings["frictional_contact_displacement_relative_tolerance"])
        conv_params.AddValue("frictional_contact_displacement_absolute_tolerance", contact_settings["frictional_contact_displacement_absolute_tolerance"])
        conv_params.AddValue("frictional_contact_residual_relative_tolerance", contact_settings["frictional_contact_residual_relative_tolerance"])
        conv_params.AddValue("frictional_contact_residual_absolute_tolerance", contact_settings["frictional_contact_residual_absolute_tolerance"])
        conv_params.AddValue("mortar_type", contact_settings["mortar_type"])
        conv_params.AddValue("condn_convergence_criterion", contact_settings["condn_convergence_criterion"])
        conv_params.AddValue("print_convergence_criterion", contact_settings["print_convergence_criterion"])
        conv_params.AddValue("ensure_contact", contact_settings["ensure_contact"])
        conv_params.AddValue("frictional_decomposed", contact_settings["frictional_decomposed"])
        conv_params.AddValue("compute_dynamic_factor", contact_settings["compute_dynamic_factor"])
        conv_params.AddValue("gidio_debug", contact_settings["gidio_debug"])

        return conv_params


def AuxiliarCreateLinearSolver(main_model_part, settings, contact_settings, linear_solver_settings, linear_solver):
    if contact_settings["rescale_linear_solver"].GetBool():
        linear_solver = KM.ScalingSolver(linear_solver, False)
    mortar_type = contact_settings["mortar_type"].GetString()
    if mortar_type == "ALMContactFrictional" or mortar_type == "ALMContactFrictionlessComponents":
        if contact_settings["use_mixed_ulm_solver"].GetBool():
            print_on_rank_zero("::[Contact Mechanical Solver]:: ", "Using MixedULMLinearSolver, definition of ALM parameters recommended")
            name_mixed_solver = contact_settings["mixed_ulm_solver_parameters"]["solver_type"].GetString()
            if name_mixed_solver == "mixed_ulm_linear_solver":
                linear_solver_name = settings["linear_solver_settings"]["solver_type"].GetString()
                if linear_solver_name == "amgcl" or linear_solver_name == "AMGCL" or linear_solver_name == "AMGCLSolver":
                    amgcl_param = KM.Parameters("""
                    {
                        "solver_type"                    : "amgcl",
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
                    amgcl_param["block_size"].SetInt(main_model_part.ProcessInfo[KM.DOMAIN_SIZE])
                    linear_solver_settings.RecursivelyValidateAndAssignDefaults(amgcl_param)
                    linear_solver = KM.AMGCLSolver(linear_solver_settings)
                mixed_ulm_solver = CSMA.MixedULMLinearSolver(linear_solver, contact_settings["mixed_ulm_solver_parameters"])
                return mixed_ulm_solver
            else:
                print_on_rank_zero("::[Contact Mechanical Solver]:: ", "Mixed solver not available: " + name_mixed_solver + ". Using not mixed linear solver")
                return linear_solver
        else:
            return linear_solver
    else:
        return linear_solver


def AuxiliarLineSearch(computing_model_part, mechanical_scheme, linear_solver, mechanical_convergence_criterion, builder_and_solver, settings, contact_settings, processes_list, post_process):
    newton_parameters = KM.Parameters("""{}""")
    return CSMA.LineSearchContactStrategy(computing_model_part,
                                            mechanical_scheme,
                                            linear_solver,
                                            mechanical_convergence_criterion,
                                            builder_and_solver,
                                            settings["max_iteration"].GetInt(),
                                            settings["compute_reactions"].GetBool(),
                                            settings["reform_dofs_at_each_step"].GetBool(),
                                            settings["move_mesh_flag"].GetBool(),
                                            newton_parameters
                                            )


def AuxiliarNewton(computing_model_part, mechanical_scheme, linear_solver, mechanical_convergence_criterion, builder_and_solver, settings, contact_settings, processes_list, post_process):
    newton_parameters = KM.Parameters("""{}""")
    newton_parameters.AddValue("adaptative_strategy", contact_settings["adaptative_strategy"])
    newton_parameters.AddValue("split_factor", contact_settings["split_factor"])
    newton_parameters.AddValue("max_number_splits", contact_settings["max_number_splits"])
    newton_parameters.AddValue("inner_loop_iterations", contact_settings["inner_loop_iterations"])
    return CSMA.ResidualBasedNewtonRaphsonContactStrategy(computing_model_part,
                                                            mechanical_scheme,
                                                            linear_solver,
                                                            mechanical_convergence_criterion,
                                                            builder_and_solver,
                                                            settings["max_iteration"].GetInt(),
                                                            settings["compute_reactions"].GetBool(),
                                                            settings["reform_dofs_at_each_step"].GetBool(),
                                                            settings["move_mesh_flag"].GetBool(),
                                                            newton_parameters,
                                                            processes_list,
                                                            post_process
                                                            )
