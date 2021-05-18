from KratosMultiphysics import Parameters
import KratosMultiphysics.SwimmingDEMApplication.swimming_DEM_procedures as SDP
import KratosMultiphysics.SwimmingDEMApplication.swimming_DEM_analysis as swimming_DEM_analysis
BaseAnalysis = swimming_DEM_analysis.SwimmingDEMAnalysis

class PreCalculatedFluidAnalysis(BaseAnalysis):
    def __init__(self, varying_parameters = Parameters("{}")):
        BaseAnalysis.__init__(self, varying_parameters)
        # self.fluid_loader_counter = self.GetFluidLoaderCounter()

    def SetBetaParameters(self):
        BaseAnalysis.SetBetaParameters(self)
        Add = self.project_parameters.AddEmptyValue
        Add("fluid_already_calculated").SetBool(True)
        Add("alpha").SetDouble(0.01)
        Add("load_derivatives").SetBool(False)
        Add("DEM_steps_per_fluid_load_step").SetInt(10)
        Add("store_fluid_in_single_precision").SetBool(True)
        Add("store_fluid_pressure_option").SetBool(True)
        Add("do_write_results_to_hdf5").SetBool(True)

    def PerformZeroStepInitializations(self):
        BaseAnalysis.PerformZeroStepInitializations(self)
        n_nodes = len(self.fluid_model_part.Nodes)
        self.project_parameters.AddEmptyValue("prerun_fluid_file_name").SetString('/mesh_' + str(n_nodes) + '_nodes.hdf5')
        self.SetFluidLoader()

    def SetFluidLoader(self):
        import hdf5_io_tools
        self.fluid_loader = hdf5_io_tools.FluidHDF5Loader(self.project_parameters,
                                                          self.all_model_parts.Get('FluidPart'),
                                                          self.all_model_parts.Get('SpheresPart'),
                                                          self.main_path)

    def FluidSolve(self, time = 'None', solve_system = True):
        if not self.project_parameters["custom_fluid"]["fluid_already_calculated"].GetBool():
            BaseAnalysis.FluidSolve(self, time, solve_system = solve_system)
            if self.project_parameters["do_write_results_to_hdf5"].GetBool():
                self.fluid_loader.FillFluidDataStep()
        else:
            BaseAnalysis.FluidSolve(self, time, solve_system = False)
            if not self.stationarity:
                self.fluid_loader.LoadFluid(time)

    def GetFluidLoaderCounter(self):
        return SDP.Counter(steps_in_cycle = self.project_parameters["DEM_steps_per_fluid_load_step"].GetInt(),
                           beginning_step = 1,
                           is_dead = not self.project_parameters["custom_fluid"]["fluid_already_calculated"].GetBool())

    def GetRunCode(self):
        code = ""

        if self.project_parameters["custom_fluid"]["fluid_already_calculated"].GetBool():
            return code + '_precalculated_fluid'
        else:
            return code

    def TheSimulationMustGoOn(self):
        it_must_go_on = self.time <= self.final_time

        if self.project_parameters["custom_fluid"]["fluid_already_calculated"].GetBool():
            return it_must_go_on and self.fluid_loader.CanLoadMoreSteps()
        else:
            return it_must_go_on
