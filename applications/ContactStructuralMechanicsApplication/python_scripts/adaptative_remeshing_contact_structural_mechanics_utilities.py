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
    import KratosMultiphysics.MeshingApplication as MeshingApplication
    missing_meshing_dependencies = False
    missing_application = ''
except ImportError as e:
    missing_meshing_dependencies = True
    # extract name of the missing application from the error message
    import re
    missing_application = re.search(r'''.*'KratosMultiphysics\.(.*)'.*''','{0}'.format(e)).group(1)

# Import base class
import adaptative_remeshing_structural_mechanics_utilities

class AdaptativeRemeshingContactMechanicalUtilities(adaptative_remeshing_structural_mechanics_utilities.AdaptativeRemeshingMechanicalUtilities):
    """These are common utilities for adaptative remeshing (for contact)
    """

    def __init__(self):
        new_parameters = KM.Parameters("""
        {
            "penalty_normal"                      : 1.0e4,
            "penalty_tangential"                  : 1.0e4
        }
        """)
        super(AdaptativeRemeshingContactMechanicalUtilities, self).__init__()
        self.adaptative_remesh_parameters["compute_error_settings"]["compute_error_extra_parameters"].AddValue("penalty_normal", new_parameters["penalty_normal"])
        self.adaptative_remesh_parameters["compute_error_settings"]["compute_error_extra_parameters"].AddValue("penalty_tangential", new_parameters["penalty_tangential"])

    def GetDefaultParameters(self):
        return self.adaptative_remesh_parameters

    def SetDefaultParameters(self, settings):
        self.adaptative_remesh_parameters = settings

    def GetConvergenceCriteria(self, error_criteria, conv_settings):
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
