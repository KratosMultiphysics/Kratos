import KratosMultiphysics as Kratos
from KratosMultiphysics import Parameters
import KratosMultiphysics.SwimmingDEMApplication as SDEM
import KratosMultiphysics.SwimmingDEMApplication.swimming_DEM_procedures as SDP
import KratosMultiphysics.SwimmingDEMApplication.cellular_flow.ethier_benchmark_analysis as ethier_benchmark_analysis
BaseAnalysis = ethier_benchmark_analysis.EthierBenchmarkAnalysis
import os

class PouliotBenchmark2DAnalysis(BaseAnalysis):
    def __init__(self, varying_parameters = Parameters("{}")):
        BaseAnalysis.__init__(self, varying_parameters)

    def ReadFluidModelParts(self):
        model_part_io_fluid = Kratos.ModelPartIO('benchmark2D')
        os.chdir(self.main_path)
        model_part_io_fluid.ReadModelPart(self._GetFluidAnalysis().fluid_model_part)

    def GetFieldUtility(self):
        self.flow_field = SDEM.PouliotFlowField2D()
        space_time_set = SDEM.SpaceTimeSet()
        self.field_utility = SDEM.FluidFieldUtility(space_time_set, self.flow_field, 1000.0, 1e-6)
        return self.field_utility

    def PerformFinalOperations(self, time = None):
        pass

    def FluidInitialize(self):
        self.fluid_model_part = self._GetFluidAnalysis().fluid_model_part
        self._GetFluidAnalysis().vars_man=self.vars_man
        self._GetFluidAnalysis().SetFluidSolverModule()
        self._GetFluidAnalysis().AddFluidVariables()
        self.AddExtraProcessInfoVariablesToFluid()
        self.ReadFluidModelParts()
        self._GetFluidAnalysis().SetFluidBufferSizeAndAddDofs()
        SDP.AddExtraDofs(self.fluid_model_part,
                         self._GetDEMAnalysis().spheres_model_part,
                         self._GetDEMAnalysis().cluster_model_part,
                         self._GetDEMAnalysis().dem_inlet_model_part,
                         self.vars_man)
        self._GetFluidAnalysis().SetFluidSolver()
        self._GetFluidAnalysis().fluid_solver.Initialize()
        self._GetFluidAnalysis().ActivateTurbulenceModel()

    def GetDebugInfo(self):
        return SDP.Counter(is_dead = True)

    def GetVolumeDebugTool(self):
        return None
