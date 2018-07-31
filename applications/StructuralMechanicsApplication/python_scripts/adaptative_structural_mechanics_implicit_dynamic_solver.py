from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics

# Check that applications were imported in the main script
KratosMultiphysics.CheckRegisteredApplications("StructuralMechanicsApplication")

# Import applications
import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication

try:
    # Check that applications were imported in the main script
    KratosMultiphysics.CheckRegisteredApplications("MeshingApplication")
    import KratosMultiphysics.MeshingApplication as MeshingApplication
    missing_meshing_dependencies = False
    missing_application = ''
except ImportError as e:
    missing_meshing_dependencies = True
    # extract name of the missing application from the error message
    import re
    missing_application = re.search(r'''.*'KratosMultiphysics\.(.*)'.*''','{0}'.format(e)).group(1)

# Import base class file
import structural_mechanics_implicit_dynamic_solver

def CreateSolver(model, custom_settings):
    return AdaptativeImplicitMechanicalSolver(model, custom_settings)

class AdaptativeImplicitMechanicalSolver(structural_mechanics_implicit_dynamic_solver.ImplicitMechanicalSolver):
    """The structural mechanics implicit dynamic solver. (Fot adaptative remeshing)

    This class creates the mechanical solvers for implicit dynamic analysis.
    It currently supports Newmark, Bossak and dynamic relaxation schemes.

    Public member variables:
    dynamic_settings -- settings for the implicit dynamic solvers.

    See structural_mechanics_solver.py for more information.
    """
    def __init__(self, model, custom_settings):
        # Set defaults and validate custom settings.
        adaptative_remesh_parameters = KratosMultiphysics.Parameters("""
        {
            "adaptative_remesh_settings" : {
                "error_mesh_tolerance" : 5.0e-3,
                "error_mesh_constant"  : 5.0e-3,
                "error_strategy_parameters":
                {
                    "minimal_size"                        : 0.01,
                    "maximal_size"                        : 1.0,
                    "error"                               : 0.01,
                    "penalty_normal"                      : 1.0e4,
                    "penalty_tangential"                  : 1.0e4,
                    "echo_level"                          : 0,
                    "set_number_of_elements"              : false,
                    "number_of_elements"                  : 1000,
                    "average_nodal_h"                     : false
                }
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
        super(AdaptativeImplicitMechanicalSolver, self).__init__(model, custom_settings)

        if (self.settings["reform_dofs_at_each_step"].GetBool() is False):
            self.print_on_rank_zero("Reform DoFs", "DoF must be reformed each time step. Switching to True")
            self.settings["reform_dofs_at_each_step"].SetBool(True)

        self.print_on_rank_zero("::[AdaptativeImplicitMechanicalSolver]:: ", "Construction finished")

    #### Private functions ####

    def AddVariables(self):
        super(AdaptativeImplicitMechanicalSolver, self).AddVariables()
        if (missing_meshing_dependencies is False):
            self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NODAL_H)
        self.print_on_rank_zero("::[AdaptativeImplicitMechanicalSolver]:: ", "Variables ADDED")

    def get_remeshing_process(self):
        if not hasattr(self, '_remeshing_process'):
            self._remeshing_process = self._create_remeshing_process()
        return self._remeshing_process

    def _create_remeshing_process(self):
        if (self.main_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE] == 2):
            remeshing_process = MeshingApplication.MmgProcess2D(self.main_model_part, self.adaptative_remesh_parameters["remeshing_parameters"])
        else:
            remeshing_process = MeshingApplication.MmgProcess3D(self.main_model_part, self.adaptative_remesh_parameters["remeshing_parameters"])

        return remeshing_process

    def _create_convergence_criterion(self):
        error_criteria = self.settings["convergence_criterion"].GetString()
        conv_settings = self._get_convergence_criterion_settings()
        if ("_with_adaptative_remesh" in error_criteria):
            conv_settings["convergence_criterion"].SetString(error_criteria.replace("_with_adaptative_remesh", ""))
        # If we just use the adaptative convergence criteria
        if (missing_meshing_dependencies is True):
            if ("adaptative_remesh" in error_criteria):
                raise NameError('The AdaptativeErrorCriteria can not be used without compiling the MeshingApplication')
        else:
            if (error_criteria == "adaptative_remesh_criteria"):
                adaptative_error_criteria = StructuralMechanicsApplication.ErrorMeshCriteria(self.adaptative_remesh_parameters["adaptative_remesh_settings"])
                adaptative_error_criteria.SetEchoLevel(conv_settings["echo_level"].GetInt())
                return adaptative_error_criteria

        # Regular convergence criteria
        import convergence_criteria_factory
        convergence_criterion = convergence_criteria_factory.convergence_criterion(conv_settings)

        # If we combine the regular convergence criteria with adaptative
        if (missing_meshing_dependencies is False):
            if ("_with_adaptative_remesh" in error_criteria):
                adaptative_error_criteria = StructuralMechanicsApplication.ErrorMeshCriteria(self.adaptative_remesh_parameters["adaptative_remesh_settings"])
                adaptative_error_criteria.SetEchoLevel(conv_settings["echo_level"].GetInt())
                convergence_criterion.mechanical_convergence_criterion = KratosMultiphysics.AndCriteria(convergence_criterion.mechanical_convergence_criterion, adaptative_error_criteria)
        return convergence_criterion.mechanical_convergence_criterion


