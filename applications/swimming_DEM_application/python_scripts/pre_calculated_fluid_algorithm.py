from KratosMultiphysics import *
from KratosMultiphysics.DEMApplication import *
import swimming_DEM_algorithm
import swimming_DEM_procedures as SDP
import math
BaseAlgorithm = swimming_DEM_algorithm.Algorithm
import h5py

class Algorithm(BaseAlgorithm):
    def __init__(self, pp):
        BaseAlgorithm.__init__(self, pp)
        self.fluid_loader_counter = self.GetFluidLoaderCounter()

    def SetBetaParamters(self):
        BaseAlgorithm.SetBetaParamters(self)
        self.pp.CFD_DEM.alpha = 0.01
        self.pp.CFD_DEM.fluid_already_calculated = 1
        self.pp.CFD_DEM.load_derivatives = 0
        self.pp.CFD_DEM.DEM_steps_per_fluid_load_step = 1
        self.pp.CFD_DEM.store_fluid_in_single_precision = 1
        self.pp.CFD_DEM.store_fluid_pressure = 1

    def PerformZeroStepInitializations(self):
        import hdf5_io_tools
        self.fluid_loader = hdf5_io_tools.FluidHDF5Loader(self.all_model_parts.Get('FluidPart'), self.pp, self.main_path)

    def FluidSolve(self, time = 'None'):
        if not self.pp.CFD_DEM.fluid_already_calculated:
            self.fluid_solver.Solve()
            self.fluid_loader.FillFluidDataStep()

    def PerformInitialDEMStepOperations(self, time = None):
        if self.fluid_loader_counter.Tick():
            if time > 0.07:
                self.fluid_loader.LoadFluid(time)

    def GetFluidLoaderCounter(self):
        return SDP.Counter(self.pp.CFD_DEM.DEM_steps_per_fluid_load_step, 1, self.pp.CFD_DEM.fluid_already_calculated)
