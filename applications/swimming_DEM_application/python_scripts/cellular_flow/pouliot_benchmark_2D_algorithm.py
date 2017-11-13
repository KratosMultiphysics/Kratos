from KratosMultiphysics import *
from KratosMultiphysics.SwimmingDEMApplication import *
from KratosMultiphysics.DEMApplication import *
import ethier_benchmark_algorithm
import swimming_DEM_procedures as SDP
import math
import numpy as np
BaseAlgorithm = ethier_benchmark_algorithm.Algorithm

class Algorithm(BaseAlgorithm):
    def __init__(self, varying_parameters = dict()):
        BaseAlgorithm.__init__(self, varying_parameters)
        self.pp.CFD_DEM.coupling_level_type = 0 #TODO: check if this should be here in this format or writing to Json!
        self.pp.CFD_DEM.laplacian_calculation_type = 0

    def ReadFluidModelParts(self):
        model_part_io_fluid = ModelPartIO('benchmark2D')
        os.chdir(self.main_path)
        model_part_io_fluid.ReadModelPart(self.fluid_algorithm.fluid_model_part)

    def GetFieldUtility(self):
        print('PouliotFlowField2D')
        self.flow_field = PouliotFlowField2D()
        space_time_set = SpaceTimeSet()
        self.field_utility = FluidFieldUtility(space_time_set, self.flow_field, 1000.0, 1e-6)
        return self.field_utility

    def PerformFinalOperations(self, time = None):
        pass

    def FluidInitialize(self):
        self.fluid_model_part = self.fluid_algorithm.fluid_model_part
        self.fluid_algorithm.vars_man=self.vars_man
        self.fluid_algorithm.SetFluidSolverModule()
        self.fluid_algorithm.AddFluidVariables()
        self.vars_man.AddExtraProcessInfoVariablesToFluidModelPart(self.pp, self.fluid_model_part)
        self.ReadFluidModelParts()
        self.fluid_algorithm.SetFluidBufferSizeAndAddDofs()
        SDP.AddExtraDofs(self.pp, self.fluid_model_part, self.disperse_phase_algorithm.spheres_model_part, self.disperse_phase_algorithm.cluster_model_part, self.disperse_phase_algorithm.DEM_inlet_model_part)
        self.fluid_algorithm.SetFluidSolver()
        self.fluid_algorithm.fluid_solver.Initialize()
        self.fluid_algorithm.ActivateTurbulenceModel()

    def GetDebugInfo(self):
        return SDP.Counter(is_dead = True)

    def GetVolumeDebugTool(self):
        return None
