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
        Add = self.pp.CFD_DEM.AddEmptyValue
        Add("fluid_already_calculated").SetBool(True)
        Add("alpha").SetDouble(0.01)
        Add("load_derivatives").SetBool(False)
        Add("DEM_steps_per_fluid_load_step").SetInt(10)
        Add("store_fluid_in_single_precision").SetBool(True)
        Add("store_fluid_pressure_option").SetBool(True)
        Add("do_write_results_to_hdf5").SetBool(True)

    def PerformZeroStepInitializations(self):
        BaseAlgorithm.PerformZeroStepInitializations(self)
        n_nodes = len(self.fluid_model_part.Nodes)
        self.pp.CFD_DEM.AddEmptyValue("prerun_fluid_file_name").SetString('/mesh_' + str(n_nodes) + '_nodes.hdf5')
        self.SetFluidLoader()

    def SetFluidLoader(self):
        import hdf5_io_tools
        self.fluid_loader = hdf5_io_tools.FluidHDF5Loader(self.all_model_parts.Get('FluidPart'),
                                                          self.all_model_parts.Get('SpheresPart'),
                                                          self.pp,
                                                          self.main_path)

    def FluidSolve(self, time = 'None', solve_system = True):
        if not self.pp.CFD_DEM["fluid_already_calculated"].GetBool():
            BaseAlgorithm.FluidSolve(self, time, solve_system = solve_system)
            if self.pp.CFD_DEM["do_write_results_to_hdf5"].GetBool():
                self.fluid_loader.FillFluidDataStep()
        else:
            BaseAlgorithm.FluidSolve(self, time, solve_system = False)
            if not self.stationarity:
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

    def TheSimulationMustGoOn(self):
        it_must_go_on = self.time <= self.end_time

        if self.pp.CFD_DEM["fluid_already_calculated"].GetBool():
            return it_must_go_on and self.fluid_loader.CanLoadMoreSteps()
        else:
            return it_must_go_on
