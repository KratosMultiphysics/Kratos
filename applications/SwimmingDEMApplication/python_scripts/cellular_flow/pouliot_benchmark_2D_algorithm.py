from KratosMultiphysics import *
from KratosMultiphysics.SwimmingDEMApplication import *
from KratosMultiphysics.DEMApplication import *
import swimming_DEM_procedures as SDP
import math
import numpy as np
import ethier_benchmark_algorithm
BaseAnalysis = ethier_benchmark_algorithm.EthierBenchmarkAnalysis

class PouliotBenchmark2DAnalysis(BaseAnalysis):
    def __init__(self, varying_parameters = dict()):
        BaseAnalysis.__init__(self, varying_parameters)
        self.pp.CFD_DEM.coupling_level_type = 0 #TODO: check if this should be here in this format or writing to Json!
        self.pp.CFD_DEM.laplacian_calculation_type = 0

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
        SDP.AddExtraDofs(self.pp, self.fluid_model_part, self.disperse_phase_solution.spheres_model_part, self.disperse_phase_solution.cluster_model_part, self.disperse_phase_solution.DEM_inlet_model_part)
        self.fluid_solution.SetFluidSolver()
        self.fluid_solution.fluid_solver.Initialize()
        self.fluid_solution.ActivateTurbulenceModel()

    def GetDebugInfo(self):
        return SDP.Counter(is_dead = True)

    def GetVolumeDebugTool(self):
        return None
