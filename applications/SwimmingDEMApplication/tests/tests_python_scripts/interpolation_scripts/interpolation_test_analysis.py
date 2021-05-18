from KratosMultiphysics import Parameters
import KratosMultiphysics.SwimmingDEMApplication.swimming_DEM_procedures as SDP
import os
file_path = os.path.abspath(__file__)
dir_path = os.path.dirname(file_path)
from KratosMultiphysics.SwimmingDEMApplication.swimming_DEM_analysis import SwimmingDEMAnalysis

class InterpolationTestAnalysis(SwimmingDEMAnalysis):
    def __init__(self, model, varying_parameters = Parameters("{}")):
        super(InterpolationTestAnalysis, self).__init__(model, varying_parameters)
        self._GetDEMAnalysis().mdpas_folder_path = os.path.join(self._GetDEMAnalysis().main_path, 'interpolation_tests')

    def GetDebugInfo(self):
        return SDP.Counter(is_dead = 1)

    def _CreateSolver(self):
        import tests_python_scripts.interpolation_scripts.interpolation_test_solver as sdem_solver
        return sdem_solver.InterpolationTestSolver(self.model,
                                                   self.project_parameters,
                                                   self.GetFieldUtility(),
                                                   self._GetFluidAnalysis()._GetSolver(),
                                                   self._GetDEMAnalysis()._GetSolver(),
                                                   self.vars_man)

    def FinalizeSolutionStep(self):
        super().FinalizeSolutionStep()

