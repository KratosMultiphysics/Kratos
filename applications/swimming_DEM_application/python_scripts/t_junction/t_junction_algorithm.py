from KratosMultiphysics import *
from KratosMultiphysics.DEMApplication import *
import pre_calculated_fluid_algorithm
BaseAlgorithm = pre_calculated_fluid_algorithm.Algorithm
import h5py

class Algorithm(BaseAlgorithm):
    def __init__(self, varying_parameters = Parameters("{}")):
        BaseAlgorithm.__init__(self, varying_parameters)

    def PerformFinalOperations(self, time = None):
        BaseAlgorithm.PerformFinalOperations(self, time)
        self.fluid_loader.RecordParticlesInBox(self.disperse_phase_algorithm.spheres_model_part)
