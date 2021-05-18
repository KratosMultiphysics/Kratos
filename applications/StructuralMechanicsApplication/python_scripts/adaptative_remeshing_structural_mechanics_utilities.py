# Importing the Kratos Library
import KratosMultiphysics

# Import applications
import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication
from KratosMultiphysics.StructuralMechanicsApplication import convergence_criteria_factory

try:
    import KratosMultiphysics.MeshingApplication as MeshingApplication
    missing_meshing_dependencies = False
    missing_application = ''
except ImportError as e:
    missing_meshing_dependencies = True
    # extract name of the missing application from the error message
    import re
    missing_application = re.search(r'''.*'KratosMultiphysics\.(.*)'.*''','{0}'.format(e)).group(1)

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
                "framework"                            : "Lagrangian",
                "internal_variables_parameters"        :
                {
                    "allocation_size"                      : 1000,
                    "bucket_size"                          : 4,
                    "search_factor"                        : 2,
                    "interpolation_type"                   : "LST",
                    "internal_variable_interpolation_list" :[]
                },
                "save_external_files"              : false,
                "max_number_of_searchs"            : 1000,
                "echo_level"                       : 0
            }
        }
        """)

    def GetDefaultParameters(self):
        return self.adaptative_remesh_parameters

    def SetDefaultParameters(self, settings):
        self.adaptative_remesh_parameters = settings

    def GetConvergenceCriteria(self, error_criteria, conv_settings):
        if ("_with_adaptative_remesh" in error_criteria):
            conv_settings["convergence_criterion"].SetString(error_criteria.replace("_with_adaptative_remesh", ""))
        # If we just use the adaptative convergence criteria
        if (missing_meshing_dependencies is True):
            if ("adaptative_remesh" in error_criteria):
                raise NameError('The AdaptativeErrorCriteria can not be used without compiling the MeshingApplication')
        else:
            if (error_criteria == "adaptative_remesh_criteria"):
                adaptative_error_criteria = StructuralMechanicsApplication.ErrorMeshCriteria(self.adaptative_remesh_parameters["compute_error_settings"])
                adaptative_error_criteria.SetEchoLevel(conv_settings["echo_level"].GetInt())
                return adaptative_error_criteria

        # Regular convergence criteria
        convergence_criterion = convergence_criteria_factory.convergence_criterion(conv_settings)

        # If we combine the regular convergence criteria with adaptative
        if (missing_meshing_dependencies is False):
            if ("_with_adaptative_remesh" in error_criteria):
                adaptative_error_criteria = StructuralMechanicsApplication.ErrorMeshCriteria(self.adaptative_remesh_parameters["compute_error_settings"])
                adaptative_error_criteria.SetEchoLevel(conv_settings["echo_level"].GetInt())
                convergence_criterion.mechanical_convergence_criterion = KratosMultiphysics.AndCriteria(convergence_criterion.mechanical_convergence_criterion, adaptative_error_criteria)
        return convergence_criterion.mechanical_convergence_criterion
