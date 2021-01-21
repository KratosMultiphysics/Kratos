# Importing the Kratos Library
import KratosMultiphysics as KM

# Import applications
import KratosMultiphysics.ContactStructuralMechanicsApplication as CSMA

# Import base class
import KratosMultiphysics.StructuralMechanicsApplication.adaptive_remeshing.adaptative_remeshing_structural_mechanics_utilities as adaptative_remeshing_structural_mechanics_utilities

# Import factories
from KratosMultiphysics.ContactStructuralMechanicsApplication.contact_convergence_criteria_factory import ContactConvergenceCriteriaFactory

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

    def GetConvergenceCriteria(self, main_model_part, error_criteria, conv_settings, compute_error_settings):
        if "_with_adaptative_remesh" in error_criteria:
            conv_settings["convergence_criterion"].SetString(error_criteria.replace("_with_adaptative_remesh", ""))
        convergence_criterion_created = ContactConvergenceCriteriaFactory(main_model_part, conv_settings)

        # If we just use the adaptative convergence criteria
        if error_criteria == "adaptative_remesh_criteria":
            adaptative_error_criteria = CSMA.ContactErrorMeshCriteria(compute_error_settings)
            convergence_criterion_created.mechanical_convergence_criterion = KM.AndCriteria(convergence_criterion_created.GetMortarCriteria(False), adaptative_error_criteria)
        elif "with_adaptative_remesh" in error_criteria: # If we combine the regular convergence criteria with adaptative
            adaptative_error_criteria = CSMA.ContactErrorMeshCriteria(compute_error_settings)
            convergence_criterion_created.mechanical_convergence_criterion = KM.AndCriteria(convergence_criterion_created.mechanical_convergence_criterion, adaptative_error_criteria)

        return convergence_criterion_created.mechanical_convergence_criterion
