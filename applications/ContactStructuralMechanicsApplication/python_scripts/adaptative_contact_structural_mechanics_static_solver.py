from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics as KM

# Check that applications were imported in the main script
KM.CheckRegisteredApplications("StructuralMechanicsApplication")
KM.CheckRegisteredApplications("ContactStructuralMechanicsApplication")

# Import applications
import KratosMultiphysics.StructuralMechanicsApplication as SMA
import KratosMultiphysics.ContactStructuralMechanicsApplication as CSMA

try:
    # Check that applications were imported in the main script
    KM.CheckRegisteredApplications("MeshingApplication")
    import KratosMultiphysics.MeshingApplication as MA
    missing_meshing_dependencies = False
    missing_application = ''
except ImportError as e:
    missing_meshing_dependencies = True
    # extract name of the missing application from the error message
    import re
    missing_application = re.search(r'''.*'KratosMultiphysics\.(.*)'.*''','{0}'.format(e)).group(1)

# Import base class file
import contact_structural_mechanics_static_solver

def CreateSolver(model, custom_settings):
    return AdaptativeContactStaticMechanicalSolver(model, custom_settings)

class AdaptativeContactStaticMechanicalSolver(contact_structural_mechanics_static_solver.ContactStaticMechanicalSolver):
    """The structural mechanics static solver. (Fot adaptative remeshing)

    This class creates the mechanical solvers for static analysis. It currently
    supports line search, linear, arc-length, form-finding and Newton-Raphson
    strategies.

    Public member variables:
    arc_length_settings -- settings for the arc length method.

    See structural_mechanics_solver.py for more information.
    """
    def __init__(self, model, custom_settings):
        # Set defaults and validate custom settings.
        adaptative_remesh_parameters = KM.Parameters("""
        {
            "compute_error_settings" : {
                "error_mesh_tolerance" : 5.0e-3,
                "error_mesh_constant"  : 5.0e-3,
                "compute_error_extra_parameters":
                {
                    "penalty_normal"                      : 1.0e4,
                    "penalty_tangential"                  : 1.0e4,
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

        # Validate the remaining settings in the base class.
        self.validate_and_transfer_matching_settings(custom_settings, adaptative_remesh_parameters)
        self.adaptative_remesh_parameters = adaptative_remesh_parameters

        # Construct the base solver.
        super(AdaptativeContactStaticMechanicalSolver, self).__init__(model, custom_settings)
        self.print_on_rank_zero("::[AdaptativeContactStaticMechanicalSolver]:: ", "Construction finished")

    #### Private functions ####

    def AddVariables(self):
        super(AdaptativeContactStaticMechanicalSolver, self).AddVariables()
        if (missing_meshing_dependencies is False):
            self.main_model_part.AddNodalSolutionStepVariable(KM.NODAL_H)
        self.print_on_rank_zero("::[AdaptativeContactStaticMechanicalSolver]:: ", "Variables ADDED")

    def get_remeshing_process(self):
        if not hasattr(self, '_remeshing_process'):
            self._remeshing_process = self._create_remeshing_process()
        return self._remeshing_process

    def _create_remeshing_process(self):
        if (self.main_model_part.ProcessInfo[KM.DOMAIN_SIZE] == 2):
            remeshing_process = MA.MmgProcess2D(self.main_model_part, self.adaptative_remesh_parameters["remeshing_parameters"])
        else:
            remeshing_process = MA.MmgProcess3D(self.main_model_part, self.adaptative_remesh_parameters["remeshing_parameters"])

        return remeshing_process

    def get_metric_process(self):
        if not hasattr(self, '_metric_process'):
            self._metric_process = self._create_metric_process()
        return self._metric_process

    def _create_metric_process(self):
        if (self.main_model_part.ProcessInfo[KM.DOMAIN_SIZE] == 2):
            metric_process = MA.MetricErrorProcess2D(self.main_model_part, self.adaptative_remesh_parameters["metric_error_parameters"])
        else:
            metric_process = MA.MetricErrorProcess3D(self.main_model_part, self.adaptative_remesh_parameters["metric_error_parameters"])

        return metric_process

    def _create_convergence_criterion(self):
        error_criteria = self.settings["convergence_criterion"].GetString()
        conv_settings = self._get_convergence_criterion_settings()
        if ("_with_adaptative_remesh" in error_criteria):
            conv_settings["convergence_criterion"].SetString(error_criteria.replace("_with_adaptative_remesh", ""))
        import contact_convergence_criteria_factory
        convergence_criterion = contact_convergence_criteria_factory.convergence_criterion(conv_settings)

        # If we just use the adaptative convergence criteria
        if (missing_meshing_dependencies is True):
            if ("adaptative_remesh" in error_criteria):
                raise NameError('The AdaptativeErrorCriteria can not be used without compiling the MeshingApplication')
        else:
            if (error_criteria == "adaptative_remesh_criteria"):
                adaptative_error_criteria = CSMA.ContactErrorMeshCriteria(self.adaptative_remesh_parameters["compute_error_settings"])
                convergence_criterion.mechanical_convergence_criterion = KM.AndCriteria(convergence_criterion.GetMortarCriteria(False), adaptative_error_criteria)
            elif ("with_adaptative_remesh" in error_criteria): # If we combine the regular convergence criteria with adaptative
                adaptative_error_criteria = CSMA.ContactErrorMeshCriteria(self.adaptative_remesh_parameters["compute_error_settings"])
                convergence_criterion.mechanical_convergence_criterion = KM.AndCriteria(convergence_criterion.mechanical_convergence_criterion, adaptative_error_criteria)

        return convergence_criterion.mechanical_convergence_criterion

        # If we combine the regular convergence criteria with adaptative
        if (missing_meshing_dependencies is False):
            if ("with_adaptative_remesh" in error_criteria):
                adaptative_error_criteria = CSMA.ContactErrorMeshCriteria(self.adaptative_remesh_parameters["compute_error_settings"])
                convergence_criterion.mechanical_convergence_criterion = KM.AndCriteria(convergence_criterion.mechanical_convergence_criterion, adaptative_error_criteria)
        return convergence_criterion.mechanical_convergence_criterion


