from KratosMultiphysics import Model, Parameters
import KratosMultiphysics.SwimmingDEMApplication.swimming_DEM_procedures as SDP
import os
file_path = os.path.abspath(__file__)
dir_path = os.path.dirname(file_path)
from KratosMultiphysics.SwimmingDEMApplication.swimming_DEM_analysis import SwimmingDEMAnalysis

class BackwardCouplingTestAnalysis(SwimmingDEMAnalysis):
    def __init__(self, model, varying_parameters = Parameters("{}")):
        super(BackwardCouplingTestAnalysis, self).__init__(model, varying_parameters)
        self._GetDEMAnalysis().mdpas_folder_path = os.path.join(self._GetDEMAnalysis().main_path, 'backward_coupling_tests/')

    def GetDebugInfo(self):
        return SDP.Counter(is_dead = 1)


    def FinalizeSolutionStep(self):
        super(BackwardCouplingTestAnalysis, self).FinalizeSolutionStep()
