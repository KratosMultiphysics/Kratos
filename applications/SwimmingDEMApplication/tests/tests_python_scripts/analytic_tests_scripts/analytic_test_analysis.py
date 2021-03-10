import KratosMultiphysics as Kratos
from KratosMultiphysics import Parameters
import KratosMultiphysics.SwimmingDEMApplication.swimming_DEM_procedures as SDP
import os
import math
file_path = os.path.abspath(__file__)
dir_path = os.path.dirname(file_path)

from KratosMultiphysics.SwimmingDEMApplication.swimming_DEM_analysis import SwimmingDEMAnalysis

class MultipleGhostsTestAnalysis(SwimmingDEMAnalysis):
    def __init__(self, model, varying_parameters = Parameters("{}")):
        super().__init__(model, varying_parameters)

    def FinalizeSolutionStep(self):
        super().FinalizeSolutionStep()
        self.CheckTotalNumberOfCrossingParticles(self.time)

    def CheckTotalNumberOfCrossingParticles(self, time):
        import h5py

        if time > 0.0002:
            with h5py.File('flux_data_new.hdf5', 'a') as f:
                n_accum_h5 = f['/'+ str(2) +'/n_accum'][-1]
                assert n_accum_h5 == -3