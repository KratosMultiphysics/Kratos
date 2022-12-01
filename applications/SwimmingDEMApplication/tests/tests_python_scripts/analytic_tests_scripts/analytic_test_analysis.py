from KratosMultiphysics import Parameters
import os
file_path = os.path.abspath(__file__)
dir_path = os.path.dirname(file_path)
print(dir_path)
import KratosMultiphysics.KratosUnittest as KratosUnittest
from KratosMultiphysics.SwimmingDEMApplication.swimming_DEM_analysis import SwimmingDEMAnalysis

class MultipleGhostsTestAnalysis(SwimmingDEMAnalysis, KratosUnittest.TestCase):
    def __init__(self, model, varying_parameters = Parameters("{}")):
        super().__init__(model, varying_parameters)

    def FinalizeSolutionStep(self):
        super().FinalizeSolutionStep()
        self.CheckTotalNumberOfCrossingParticles()

    def CheckTotalNumberOfCrossingParticles(self):
        import h5py

        if self.time > 0.0002:
            with h5py.File('flux_data_new.hdf5', 'a') as f:
                n_accum_h5 = f['/'+ str(2) +'/n_accum'][-1]
                self.assertEqual(n_accum_h5,-3)
