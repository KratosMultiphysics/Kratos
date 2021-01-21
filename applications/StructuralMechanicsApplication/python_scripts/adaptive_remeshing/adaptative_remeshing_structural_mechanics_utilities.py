from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics

# Import applications
import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication

# Import convergence criteria factory
from KratosMultiphysics.StructuralMechanicsApplication.convergence_criteria_factory import convergence_criterion

class AdaptativeRemeshingMechanicalUtilities(object):
    """These are common utilities for adaptative remeshing
    """

    def __init__(self):
        self.adaptative_remesh_parameters = KratosMultiphysics.Parameters("""
        {
            "compute_error_settings" : {
                "error_mesh_tolerance" : 5.0e-3,
                "error_mesh_constant"  : 5.0e-3,
                "compute_error_extra_parameters":
                {
                    "stress_vector_variable"              : "CAUCHY_STRESS_VECTOR",
                    "echo_level"                          : 0
                }
            },
            "metric_error_parameters" :
            {
                "minimal_size"                        : 0.01,
                "maximal_size"                        : 1.0,
                "error_strategy_parameters":
                {
                    "target_error"                        : 0.01,
                    "set_target_number_of_elements"       : false,
                    "target_number_of_elements"           : 1000,
                    "perform_nodal_h_averaging"           : false
                },
                "echo_level"                          : 0
            },
            "remeshing_parameters":
            {
                "filename"                             : "out",
                "discretization_type"                  : "Standard",
                "isosurface_parameters"                :
                {
                    "isosurface_variable"              : "DISTANCE",
                    "nonhistorical_variable"           : false,
                    "remove_internal_regions"          : false
                },
                "framework"                            : "Lagrangian",
                "internal_variables_parameters"        :
                {
                    "allocation_size"                      : 1000,
                    "bucket_size"                          : 4,
                    "search_factor"                        : 2,
                    "interpolation_type"                   : "LST",
                    "internal_variable_interpolation_list" :[]
                },
                "force_sizes"                          :
                {
                    "force_min"                           : false,
                    "minimal_size"                        : 0.1,
                    "force_max"                           : false,
                    "maximal_size"                        : 10.0
                },
                "advanced_parameters"                  :
                {
                    "force_hausdorff_value"               : false,
                    "hausdorff_value"                     : 0.0001,
                    "no_move_mesh"                        : false,
                    "no_surf_mesh"                        : false,
                    "no_insert_mesh"                      : false,
                    "no_swap_mesh"                        : false,
                    "normal_regularization_mesh"          : false,
                    "deactivate_detect_angle"             : false,
                    "force_gradation_value"               : false,
                    "gradation_value"                     : 1.3
                },
                "collapse_prisms_elements"             : false,
                "save_external_files"                  : false,
                "save_colors_files"                    : false,
                "save_mdpa_file"                       : false,
                "max_number_of_searchs"                : 1000,
                "preserve_flags"                       : true,
                "interpolate_non_historical"           : true,
                "extrapolate_contour_values"           : true,
                "surface_elements"                     : false,
                "search_parameters"                    : {
                    "allocation_size"                     : 1000,
                    "bucket_size"                         : 4,
                    "search_factor"                       : 2.0
                },
                "echo_level"                           : 0,
                "debug_result_mesh"                    : false,
                "step_data_size"                       : 0,
                "initialize_entities"                  : true,
                "remesh_at_non_linear_iteration"       : false,
                "buffer_size"                          : 0

            }
        }
        """)

    def GetDefaultParameters(self):
        return self.adaptative_remesh_parameters

    def SetDefaultParameters(self, settings):
        self.adaptative_remesh_parameters = settings

    def GetConvergenceCriteria(self, error_criteria, conv_settings, compute_error_settings):
        if "_with_adaptative_remesh" in error_criteria:
            conv_settings["convergence_criterion"].SetString(error_criteria.replace("_with_adaptative_remesh", ""))
        # If we just use the adaptative convergence criteria
        if error_criteria == "adaptative_remesh_criteria":
            adaptative_error_criteria = StructuralMechanicsApplication.ErrorMeshCriteria(compute_error_settings)
            adaptative_error_criteria.SetEchoLevel(conv_settings["echo_level"].GetInt())
            return adaptative_error_criteria

        # Regular convergence criteria
        convergence_criterion_created = convergence_criterion(conv_settings)

        # If we combine the regular convergence criteria with adaptative
        if "_with_adaptative_remesh" in error_criteria:
            adaptative_error_criteria = StructuralMechanicsApplication.ErrorMeshCriteria(compute_error_settings)
            adaptative_error_criteria.SetEchoLevel(conv_settings["echo_level"].GetInt())
            convergence_criterion_created.mechanical_convergence_criterion = KratosMultiphysics.AndCriteria(convergence_criterion_created.mechanical_convergence_criterion, adaptative_error_criteria)
        return convergence_criterion_created.mechanical_convergence_criterion
