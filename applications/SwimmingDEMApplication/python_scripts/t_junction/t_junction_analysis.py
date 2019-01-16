from KratosMultiphysics import *
from KratosMultiphysics.DEMApplication import *
from KratosMultiphysics.SwimmingDEMApplication import *
import pre_calculated_fluid_analysis
BaseAnalysis = pre_calculated_fluid_analysis.PreCalculatedFluidAnalysis

class TJunctionAnalysis(BaseAnalysis):
    def __init__(self, varying_parameters = Parameters("{}")):
        BaseAnalysis.__init__(self, varying_parameters)
        final_time = self.pp.CFD_DEM.AddEmptyValue("FinalTime").GetDouble()
        L = 0.0048 # the channel width
        center_x = 0.0044
        self.bbox_watcher = BoundingBoxRule(0.0, 2 * final_time,
                                            center_x - L, center_x + L,
                                            -0.007, -0.002,
                                            -0.005, 0.001)

    def PerformFinalOperations(self, time = None):
        self.particles_loader.RecordParticlesInBox(self.bbox_watcher)
        BaseAnalysis.PerformFinalOperations(self, time)

    def PerformZeroStepInitializations(self):
        BaseAnalysis.PerformZeroStepInitializations(self)
        import hdf5_io_tools
        self.particles_loader = hdf5_io_tools.ParticleHistoryLoader(self.all_model_parts.Get('SpheresPart'), self.disperse_phase_solution.watcher, self.pp, self.main_path)

    def FluidSolve(self, time = 'None', solve_system = True):
        BaseAnalysis.FluidSolve(self, time, solve_system)
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
