# importing the Kratos Library
import KratosMultiphysics as KM
import KratosMultiphysics.ShallowWaterApplication as SW

## Import base class file
from KratosMultiphysics.ShallowWaterApplication.shallow_water_base_solver import ShallowWaterBaseSolver
from KratosMultiphysics.ShallowWaterApplication.wave_solver import WaveSolver

def CreateSolver(model, custom_settings):
    return LagrangianShallowWaterSolver(model, custom_settings)

class LagrangianShallowWaterSolver(ShallowWaterBaseSolver):
    def __init__(self, model, settings):
        super().__init__(model, settings)
        self.min_buffer_size = 1
        self.mesh_solver = WaveSolver(self.model, self.settings["mesh_solver_settings"])

        eulerian_model_part = self.main_model_part
        lagrangian_model_part = self.mesh_solver.main_model_part
        self.mesh_moving = SW.MoveMeshUtility(lagrangian_model_part, eulerian_model_part, self.settings["mesh_moving_settings"])

    def AddVariables(self):
        super().AddVariables()
        self.main_model_part.AddNodalSolutionStepVariable(KM.ACCELERATION)
        self.main_model_part.AddNodalSolutionStepVariable(SW.VERTICAL_VELOCITY)
        self.mesh_solver.AddVariables()

    def ImportModelPart(self):
        pass

    def PrepareModelPart(self):
        self.mesh_solver.PrepareModelPart()

    def AddDofs(self):
        self.mesh_solver.AddDofs()

    def AdvanceInTime(self, current_time):
        return self.mesh_solver.AdvanceInTime(current_time)

    def GetComputingModelPart(self):
        return self.mesh_solver.GetComputingModelPart()

    def Initialize(self):
        self.mesh_moving.Initialize()
        self.mesh_solver.Initialize()

    def InitializeSolutionStep(self):
        self.mesh_moving.MoveMesh()
        self.mesh_solver.InitializeSolutionStep()

    def Predict(self):
        self.mesh_solver.Predict()

    def SolveSolutionStep(self):
        self.mesh_solver.SolveSolutionStep()

    def FinalizeSolutionStep(self):
        self.mesh_solver.FinalizeSolutionStep()
        self.mesh_moving.MapResults()

    def Finalize(self):
        self.mesh_solver.Finalize()

    def Check(self):
        self.mesh_solver.Check()
        self.mesh_moving.Check()

    def Clear(self):
        self.mesh_solver.Clear()

    @classmethod
    def GetDefaultParameters(cls):
        default_settings = KM.Parameters("""
        {
            "solver_type"              : "shallow_water_base_solver",
            "model_part_name"          : "eulerian_model_part",
            "domain_size"              : 2,
            "echo_level"               : 0,
            "mesh_solver_settings"     : {},
            "mesh_moving_settings"     : {}
        }""")
        return default_settings
