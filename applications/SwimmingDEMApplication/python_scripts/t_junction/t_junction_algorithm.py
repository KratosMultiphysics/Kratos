from KratosMultiphysics import *
from KratosMultiphysics.DEMApplication import *
from KratosMultiphysics.SwimmingDEMApplication import *
import pre_calculated_fluid_algorithm
BaseAlgorithm = pre_calculated_fluid_algorithm.Algorithm
import h5py

class Algorithm(BaseAlgorithm):
    def __init__(self, varying_parameters = Parameters("{}")):
        BaseAlgorithm.__init__(self, varying_parameters)
        end_time = self.pp.CFD_DEM.AddEmptyValue("FinalTime").GetDouble()
        L = 0.0048 # the channel width
        center_x = 0.0044
        self.bbox_watcher = BoundingBoxRule(0.0, 2 * end_time,
                                            center_x - L, center_x + L,
                                            -0.007, -0.002,
                                            -0.005, 0.001)

    def PerformFinalOperations(self, time = None):
        self.particles_loader.RecordParticlesInBox(self.bbox_watcher)
        BaseAlgorithm.PerformFinalOperations(self, time)

    def PerformZeroStepInitializations(self):
        BaseAlgorithm.PerformZeroStepInitializations(self)
        import hdf5_io_tools
        self.particles_loader = hdf5_io_tools.ParticleHistoryLoader(self.all_model_parts.Get('SpheresPart'), self.disperse_phase_solution.watcher, self.pp, self.main_path)

    def FluidSolve(self, time = 'None', solve_system = True):
        BaseAlgorithm.FluidSolve(self, time, solve_system)
        self.particles_loader.UpdateListOfAllParticles()

    def ModifyResultsFolderName(self, time):
        import os
        new_path = self.post_path + '_' + self.particles_loader.GetRunCode()
        if os.path.exists(new_path):
            import shutil
            shutil.rmtree(new_path)
        os.rename(self.post_path, new_path)

    def GetRunCode(self):
        return ''
