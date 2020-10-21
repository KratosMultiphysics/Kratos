from KratosMultiphysics import Model, Parameters
import KratosMultiphysics.SwimmingDEMApplication.swimming_DEM_procedures as SDP
from KratosMultiphysics.SwimmingDEMApplication.swimming_DEM_analysis import SwimmingDEMAnalysis
import os

class BackwardCouplingTestAnalysis(SwimmingDEMAnalysis):
    def __init__(self, model, varying_parameters = Parameters("{}")):
        super(BackwardCouplingTestAnalysis, self).__init__(model, varying_parameters)
        self._GetDEMAnalysis().mdpas_folder_path = os.path.join(self._GetDEMAnalysis().main_path, 'interpolation_tests')

    def GetDebugInfo(self):
        return SDP.Counter(is_dead = 1)

    def FinalizeSolutionStep(self):
        super(BackwardCouplingTestAnalysis, self).FinalizeSolutionStep()
