from KratosMultiphysics import *
from KratosMultiphysics.SwimmingDEMApplication import *
from KratosMultiphysics.DEMApplication import *
import swimming_DEM_procedures as SDP
import ethier_benchmark_analysis
BaseAnalysis = ethier_benchmark_analysis.EthierBenchmarkAnalysis

class PouliotBenchmark2DAnalysis(BaseAnalysis):
    def __init__(self, varying_parameters = Parameters("{}")):
        BaseAnalysis.__init__(self, varying_parameters)

    def ReadFluidModelParts(self):
        model_part_io_fluid = ModelPartIO('benchmark2D')
        os.chdir(self.main_path)
        model_part_io_fluid.ReadModelPart(self.fluid_solution.fluid_model_part)

    def GetFieldUtility(self):
        print('PouliotFlowField2D')
        self.flow_field = PouliotFlowField2D()
        space_time_set = SpaceTimeSet()
        self.field_utility = FluidFieldUtility(space_time_set, self.flow_field, 1000.0, 1e-6)
        return self.field_utility

    def PerformFinalOperations(self, time = None):
        pass

    def FluidInitialize(self):
        self.fluid_model_part = self.fluid_solution.fluid_model_part
        self.fluid_solution.vars_man=self.vars_man
        self.fluid_solution.SetFluidSolverModule()
        self.fluid_solution.AddFluidVariables()
        self.AddExtraProcessInfoVariablesToFluid()
        self.ReadFluidModelParts()
        self.fluid_solution.SetFluidBufferSizeAndAddDofs()
        SDP.AddExtraDofs(self.fluid_model_part,
                         self.disperse_phase_solution.spheres_model_part,
                         self.disperse_phase_solution.cluster_model_part,
                         self.disperse_phase_solution.dem_inlet_model_part,
                         self.vars_man)
        self.fluid_solution.SetFluidSolver()
        self.fluid_solution.fluid_solver.Initialize()
        self.fluid_solution.ActivateTurbulenceModel()

    def GetDebugInfo(self):
        return SDP.Counter(is_dead = True)

    def GetVolumeDebugTool(self):
        return None
