# Importing the Kratos Library
import KratosMultiphysics

# Import base class file
from KratosMultiphysics.StructuralMechanicsApplication import structural_mechanics_static_solver

# Import auxiliar methods
from KratosMultiphysics.StructuralMechanicsApplication.adaptive_remeshing import auxiliar_methods_adaptative_solvers
from KratosMultiphysics.StructuralMechanicsApplication.adaptive_remeshing.adaptative_remeshing_structural_mechanics_utilities import AdaptativeRemeshingMechanicalUtilities

def CreateSolver(model, custom_settings):
    return AdaptativeRemeshingStaticMechanicalSolver(model, custom_settings)

class AdaptativeRemeshingStaticMechanicalSolver(structural_mechanics_static_solver.StaticMechanicalSolver):
    """The structural mechanics static solver. (Fot adaptative remeshing)
    See structural_mechanics_static_solver.py for more information.
    """
    def __init__(self, model, custom_settings):
        # Construct the base solver.
        super().__init__(model, custom_settings)
        KratosMultiphysics.Logger.PrintInfo("::[AdaptativeRemeshingStaticMechanicalSolver]:: ", "Construction finished")

    #### Private functions ####

    def AddVariables(self):
        super(AdaptativeRemeshingStaticMechanicalSolver, self).AddVariables()
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NODAL_H)
        KratosMultiphysics.Logger.PrintInfo("::[AdaptativeRemeshingStaticMechanicalSolver]:: ", "Variables ADDED")

    def get_remeshing_process(self):
        if not hasattr(self, '_remeshing_process'):
            self._remeshing_process = self._create_remeshing_process()
        return self._remeshing_process

    def _create_remeshing_process(self):
        return auxiliar_methods_adaptative_solvers.CreateRemeshingProcess(self.main_model_part, self.settings)

    def get_metric_process(self):
        if not hasattr(self, '_metric_process'):
            self._metric_process = self._create_metric_process()
        return self._metric_process

    def _create_metric_process(self):
        return auxiliar_methods_adaptative_solvers.CreateMetricProcess(self.main_model_part, self.settings)

    def _create_convergence_criterion(self):
        error_criteria = self.settings["convergence_criterion"].GetString()
        conv_settings = self._get_convergence_criterion_settings()
        return AdaptativeRemeshingMechanicalUtilities().GetConvergenceCriteria(error_criteria, conv_settings, self.settings["compute_error_settings"])

    @classmethod
    def GetDefaultSettings(cls):
        # Set defaults and validate custom settings.
        this_defaults = AdaptativeRemeshingMechanicalUtilities().GetDefaultParameters()
        this_defaults.RecursivelyAddMissingParameters(super(AdaptativeRemeshingStaticMechanicalSolver, cls).GetDefaultSettings())
        return this_defaults
