# importing the Kratos Library
import KratosMultiphysics as KM
import KratosMultiphysics.ShallowWaterApplication as SW

## Import base class file
from KratosMultiphysics.python_solver import PythonSolver
from KratosMultiphysics.ShallowWaterApplication.wave_solver import WaveSolver

def CreateSolver(model, custom_settings):
    return LagrangianShallowWaterSolver(model, custom_settings)

class LagrangianShallowWaterSolver(PythonSolver):
    def __init__(self, model, settings):
        super().__init__(model, settings)

        self.mesh_solver = WaveSolver(self.model, self.settings["mesh_solver_settings"])

        # Either retrieve the model part from the model or create a new one
        model_part_name = self.settings["model_part_name"].GetString()
        if self.model.HasModelPart(model_part_name):
            self.main_model_part = self.model.GetModelPart(model_part_name)
        else:
            self.main_model_part = self.model.CreateModelPart(model_part_name)

        eulerian_model_part = self.main_model_part
        lagrangian_model_part = self.mesh_solver.main_model_part
        self.mesh_moving = SW.MoveMeshUtility(lagrangian_model_part, eulerian_model_part, self.settings["mesh_moving_settings"])

    def AddVariables(self):
        self.mesh_solver.AddVariables()
        self.mesh_solver.main_model_part.AddNodalSolutionStepVariable(KM.DISPLACEMENT)
        KM.MergeVariableListsUtility().Merge(self.mesh_solver.main_model_part, self.main_model_part)

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
        return self.mesh_solver.SolveSolutionStep()

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
            "solver_type"              : "",
            "model_part_name"          : "eulerian_model_part",
            "echo_level"               : 0,
            "mesh_solver_settings"     : {},
            "mesh_moving_settings"     : {}
        }""")
        return default_settings
