# Importing the Kratos Library
import KratosMultiphysics as KM

# Import applications
import KratosMultiphysics.StructuralMechanicsApplication as SMA
import KratosMultiphysics.ContactStructuralMechanicsApplication as CSMA

import KratosMultiphysics.kratos_utilities as kratos_utilities
if kratos_utilities.CheckIfApplicationsAvailable("MeshingApplication"):
    has_meshing_application = True
else:
    has_meshing_application = False

# Import adaptive remeshing utilities
import KratosMultiphysics.ContactStructuralMechanicsApplication.adaptative_remeshing_contact_structural_mechanics_utilities as adaptative_remeshing_contact_structural_mechanics_utilities

# Import base class file
import KratosMultiphysics.ContactStructuralMechanicsApplication.contact_structural_mechanics_implicit_dynamic_solver as contact_structural_mechanics_implicit_dynamic_solver

def CreateSolver(model, custom_settings):
    return AdaptativeRemeshingContactImplicitMechanicalSolver(model, custom_settings)

class AdaptativeRemeshingContactImplicitMechanicalSolver(contact_structural_mechanics_implicit_dynamic_solver.ContactImplicitMechanicalSolver):
    """The contact structural mechanics implicit dynamic solver. (For adaptative remeshing)
    See contact_structural_mechanics_implicit_dynamic_solver.py for more information.
    """
    def __init__(self, model, custom_settings):
        # Set defaults and validate custom settings.
        self.adaptative_remeshing_utilities = adaptative_remeshing_contact_structural_mechanics_utilities.AdaptativeRemeshingContactMechanicalUtilities()

        # Construct the base solver.
        super(AdaptativeRemeshingContactImplicitMechanicalSolver, self).__init__(model, custom_settings)
        KM.Logger.PrintInfo("::[AdaptativeRemeshingContactImplicitMechanicalSolver]:: ", "Construction finished")

    #### Private functions ####

    def AddVariables(self):
        super(AdaptativeRemeshingContactImplicitMechanicalSolver, self).AddVariables()
        if has_meshing_application:
            self.main_model_part.AddNodalSolutionStepVariable(KM.NODAL_H)
        KM.Logger.PrintInfo("::[AdaptativeRemeshingContactImplicitMechanicalSolver]:: ", "Variables ADDED")

    def get_remeshing_process(self):
        if not hasattr(self, '_remeshing_process'):
            self._remeshing_process = self._create_remeshing_process()
        return self._remeshing_process

    def _create_remeshing_process(self):
        if self.main_model_part.ProcessInfo[KM.DOMAIN_SIZE] == 2:
            remeshing_process = MA.MmgProcess2D(self.main_model_part, self.settings["remeshing_parameters"])
        else:
            remeshing_process = MA.MmgProcess3D(self.main_model_part, self.settings["remeshing_parameters"])

        return remeshing_process

    def get_metric_process(self):
        if not hasattr(self, '_metric_process'):
            self._metric_process = self._create_metric_process()
        return self._metric_process

    def _create_metric_process(self):
        if self.main_model_part.ProcessInfo[KM.DOMAIN_SIZE] == 2:
            metric_process = MA.MetricErrorProcess2D(self.main_model_part, self.settings["metric_error_parameters"])
        else:
            metric_process = MA.MetricErrorProcess3D(self.main_model_part, self.settings["metric_error_parameters"])

        return metric_process

    def _CreateConvergenceCriterion(self):
        error_criteria = self.settings["convergence_criterion"].GetString()
        conv_settings = self._get_convergence_criterion_settings()
        return self.adaptative_remeshing_utilities.GetConvergenceCriteria(error_criteria, conv_settings)

    @classmethod
    def GetDefaultParameters(cls):
        this_defaults = adaptative_remeshing_contact_structural_mechanics_utilities.AdaptativeRemeshingContactMechanicalUtilities().GetDefaultParameters()
        this_defaults.AddMissingParameters(super().GetDefaultParameters())
        return this_defaults
