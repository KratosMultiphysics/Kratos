from KratosMultiphysics import *
from KratosMultiphysics.DEMApplication import *
import swimming_DEM_algorithm
import swimming_DEM_procedures as SDP
import math
BaseAlgorithm = swimming_DEM_algorithm.Algorithm
import h5py

class Algorithm(BaseAlgorithm):
    def __init__(self, varying_parameters = Parameters("{}")):
        BaseAlgorithm.__init__(self, varying_parameters)
        # self.fluid_loader_counter = self.GetFluidLoaderCounter()

    def SetBetaParameters(self):
        BaseAlgorithm.SetBetaParameters(self)
        self.pp.CFD_DEM.AddEmptyValue("alpha").SetDouble(0.01)
        self.pp.CFD_DEM.AddEmptyValue("fluid_already_calculated").SetBool(True)
        self.pp.CFD_DEM.AddEmptyValue("load_derivatives").SetBool(False)
        self.pp.CFD_DEM.AddEmptyValue("DEM_steps_per_fluid_load_step").SetInt(10)
        self.pp.CFD_DEM.AddEmptyValue("store_fluid_in_single_precision").SetBool(True)
        self.pp.CFD_DEM.AddEmptyValue("store_fluid_pressure_option").SetBool(True)

    def PerformZeroStepInitializations(self):
        import hdf5_io_tools
        self.fluid_loader = hdf5_io_tools.FluidHDF5Loader(self.all_model_parts.Get('FluidPart'), self.pp, self.main_path)

    def FluidSolve(self, time = 'None'):
        if not self.pp.CFD_DEM["fluid_already_calculated"].GetBool():
            self.fluid_algorithm.fluid_solver.Solve()
            self.fluid_loader.FillFluidDataStep()
        else:
            self.fluid_loader.LoadFluid(time)

    def PerformInitialDEMStepOperations(self, time = None):
        pass

    def GetFluidLoaderCounter(self):
        return SDP.Counter(steps_in_cycle = self.pp.CFD_DEM["DEM_steps_per_fluid_load_step"].GetInt(),
                           beginning_step = 1,
                           is_dead = not self.pp.CFD_DEM["fluid_already_calculated"].GetBool())

    def GetRunCode(self):
        code = SDP.CreateRunCode(self.pp)

        if self.pp.CFD_DEM["fluid_already_calculated"].GetBool():
            return code + '_precalculated_fluid'
        else:
            return code

    def TheSimulationMustGoON(self):
        it_must_go_on = self.time <= self.final_time

        if self.pp.CFD_DEM["fluid_already_calculated"].GetBool():
            return it_must_go_on and self.fluid_loader.CanLoadMoreSteps()
        else:
            return it_must_go_on
