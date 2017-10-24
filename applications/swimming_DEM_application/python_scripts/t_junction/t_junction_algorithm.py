from KratosMultiphysics import *
from KratosMultiphysics.DEMApplication import *
import pre_calculated_fluid_algorithm
BaseAlgorithm = pre_calculated_fluid_algorithm.Algorithm
import h5py

class Algorithm(BaseAlgorithm):
    def __init__(self, varying_parameters = Parameters("{}")):
        BaseAlgorithm.__init__(self, varying_parameters)
        self.bb_watcher_min = [-0.012, -0.007, -0.005]
        self.bb_watcher_max = [0.0208, 0.018, 0.001]

    def PerformFinalOperations(self, time = None):
        BaseAlgorithm.PerformFinalOperations(self, time)
        self.fluid_loader.RecordParticlesInBox(particles_model_part = self.disperse_phase_algorithm.spheres_model_part,
                                               bb_low = self.bb_watcher_min,
                                               bb_high = self.bb_watcher_max)

    def PerformZeroStepInitializations(self):
        BaseAlgorithm.PerformZeroStepInitializations(self)
        import hdf5_io_tools
        self.particles_loader = hdf5_io_tools.ParticleHistoryLoader(self.all_model_parts.Get('SpheresPart'), self.pp, self.main_path)

    def FluidSolve(self, time = 'None'):
        BaseAlgorithm.FluidSolve(self, time)
        ids = []
        X0s = []
        Y0s = []
        Z0s = []
        radii = []
        times = []

        self.disperse_phase_algorithm.watcher.GetNewParticlesData(ids, X0s, Y0s, Z0s, radii, times)
        self.particles_loader.UpdateListOfAllParticles(ids, X0s, Y0s, Z0s, radii, times)
