"""Linear transient heat conduction on fixed SPH particles."""

import KratosMultiphysics as KM
import KratosMultiphysics.SPHApplication as SPH
from KratosMultiphysics.python_solver import PythonSolver
from KratosMultiphysics import python_linear_solver_factory


def CreateSolver(model, settings):
    return SPHThermalSolver(model, settings)


class SPHThermalSolver(PythonSolver):
    @classmethod
    def GetDefaultParameters(cls):
        defaults = KM.Parameters('''{
            "solver_type": "thermal",
            "model_part_name": "Thermal",
            "domain_size": 2,
            "echo_level": 0,
            "model_import_settings": {"input_type": "use_input_model_part"},
            "material_import_settings": {"materials_filename": ""},
            "time_stepping": {"time_step": 0.01},
            "transient_parameters": {"theta": 0.5},
            "linear_solver_settings": {"solver_type": "skyline_lu_factorization"}
        }''')
        defaults.AddMissingParameters(super().GetDefaultParameters())
        return defaults

    def __init__(self, model, settings):
        super().__init__(model, settings)
        name = self.settings["model_part_name"].GetString()
        self.main_model_part = (model[name] if model.HasModelPart(name)
                                else model.CreateModelPart(name))
        self.main_model_part.ProcessInfo[KM.DOMAIN_SIZE] = self.settings["domain_size"].GetInt()
        theta = self.settings["transient_parameters"]["theta"].GetDouble()
        if not 0.0 <= theta <= 1.0:
            raise ValueError("theta must be in [0, 1]")
        if self.settings["time_stepping"]["time_step"].GetDouble() <= 0.0:
            raise ValueError("time_step must be positive")
        self.main_model_part.ProcessInfo[KM.TIME_INTEGRATION_THETA] = theta
        self.main_model_part.ProcessInfo[KM.DELTA_TIME] = self.settings["time_stepping"]["time_step"].GetDouble()

    def AddVariables(self):
        self.main_model_part.AddNodalSolutionStepVariable(KM.TEMPERATURE)
        self.main_model_part.AddNodalSolutionStepVariable(KM.HEAT_FLUX)

    def ImportModelPart(self):
        self._ImportModelPart(self.main_model_part, self.settings["model_import_settings"])

    def PrepareModelPart(self):
        self.main_model_part.SetBufferSize(self.GetMinimumBufferSize())
        filename = self.settings["material_import_settings"]["materials_filename"].GetString()
        if filename:
            material_settings = KM.Parameters('{"Parameters":{"materials_filename":""}}')
            material_settings["Parameters"]["materials_filename"].SetString(filename)
            KM.ReadMaterialsUtility(material_settings, self.model)

    def GetMinimumBufferSize(self):
        return 2

    def AddDofs(self):
        for node in self.main_model_part.Nodes:
            node.AddDof(KM.TEMPERATURE)

    def GetComputingModelPart(self):
        return self.main_model_part

    def Initialize(self):
        # Volumes, neighbours and kernel corrections must already be initialized.
        linear_solver = python_linear_solver_factory.ConstructSolver(self.settings["linear_solver_settings"])
        self.strategy = KM.ResidualBasedLinearStrategy(
            self.main_model_part, KM.ResidualBasedIncrementalUpdateStaticScheme(),
            linear_solver, False, False, False, False)
        self.strategy.SetEchoLevel(self.echo_level)
        self.strategy.Initialize()

    def AdvanceInTime(self, current_time):
        new_time = current_time + self.settings["time_stepping"]["time_step"].GetDouble()
        self.main_model_part.CloneTimeStep(new_time)
        info = self.main_model_part.ProcessInfo
        info[KM.STEP] += 1
        return new_time

    def InitializeSolutionStep(self):
        self.strategy.InitializeSolutionStep()

    def Predict(self):
        self.strategy.Predict()

    def SolveSolutionStep(self):
        self.strategy.SolveSolutionStep()
        return True

    def FinalizeSolutionStep(self):
        self.strategy.FinalizeSolutionStep()

    def Check(self):
        return self.strategy.Check()

    def Clear(self):
        self.strategy.Clear()
