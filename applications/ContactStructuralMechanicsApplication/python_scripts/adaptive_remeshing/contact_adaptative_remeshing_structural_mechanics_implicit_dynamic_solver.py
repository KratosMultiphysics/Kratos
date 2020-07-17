# Importing the Kratos Library
import KratosMultiphysics as KM

import KratosMultiphysics.kratos_utilities as kratos_utilities
if kratos_utilities.CheckIfApplicationsAvailable("MeshingApplication"):
    has_meshing_application = True
    import KratosMultiphysics.MeshingApplication as MA
else:
    has_meshing_application = False

# Import base class file
from KratosMultiphysics.ContactStructuralMechanicsApplication import contact_structural_mechanics_implicit_dynamic_solver

# Import utilities
from KratosMultiphysics.ContactStructuralMechanicsApplication.adaptive_remeshing.adaptative_remeshing_contact_structural_mechanics_utilities import AdaptativeRemeshingContactMechanicalUtilities

def CreateSolver(model, custom_settings):
    return AdaptativeRemeshingContactImplicitMechanicalSolver(model, custom_settings)

class AdaptativeRemeshingContactImplicitMechanicalSolver(contact_structural_mechanics_implicit_dynamic_solver.ContactImplicitMechanicalSolver):
    """The contact structural mechanics implicit dynamic solver. (Fot adaptative remeshing)
    See contact_structural_mechanics_implicit_dynamic_solver.py for more information.
    """
    def __init__(self, model, custom_settings):
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
        if has_meshing_application:
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
        if has_meshing_application:
            if self.main_model_part.ProcessInfo[KM.DOMAIN_SIZE] == 2:
                metric_process = MA.MetricErrorProcess2D(self.main_model_part, self.settings["metric_error_parameters"])
            else:
                metric_process = MA.MetricErrorProcess3D(self.main_model_part, self.settings["metric_error_parameters"])

            return metric_process

    def _create_convergence_criterion(self):
        error_criteria = self.settings["convergence_criterion"].GetString()
        conv_settings = self._get_convergence_criterion_settings()
        return AdaptativeRemeshingContactMechanicalUtilities().GetConvergenceCriteria(self.main_model_part, error_criteria, conv_settings, self.settings["compute_error_settings"])

    @classmethod
    def GetDefaultSettings(cls):
        # Set defaults and validate custom settings.
        this_defaults = AdaptativeRemeshingContactMechanicalUtilities().GetDefaultParameters()
        this_defaults.RecursivelyAddMissingParameters(super(AdaptativeRemeshingContactImplicitMechanicalSolver, cls).GetDefaultSettings())
        return this_defaults
