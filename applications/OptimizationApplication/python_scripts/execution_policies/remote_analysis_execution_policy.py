from importlib import import_module
from typing import Any

import KratosMultiphysics as Kratos
import KratosMultiphysics.OptimizationApplication as KratosOA
from KratosMultiphysics.OptimizationApplication.execution_policies.execution_policy import ExecutionPolicy
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem import OptimizationProblem
from KratosMultiphysics.CoSimulationApplication import CoSimIO
from KratosMultiphysics.OptimizationApplication.utilities.logger_utilities import time_decorator


def Factory(model: Kratos.Model, parameters: Kratos.Parameters, optimization_problem: OptimizationProblem) -> ExecutionPolicy:
    if not parameters.Has("name"):
        raise RuntimeError(f"RemoteAnalysisExecutionPolicy instantiation requires a \"name\" in parameters [ parameters = {parameters}].")
    if not parameters.Has("settings"):
        raise RuntimeError(f"RemoteAnalysisExecutionPolicy instantiation requires a \"settings\" in parameters [ parameters = {parameters}].")

    return RemoteAnalysisExecutionPolicy(parameters["name"].GetString(), model, parameters["settings"], optimization_problem)

class RemoteAnalysisExecutionPolicy(ExecutionPolicy):
    def __init__(self, name: str, model: Kratos.Model, parameters: Kratos.Parameters, optimization_problem: OptimizationProblem):
        super().__init__(name)

        default_settings = Kratos.Parameters("""{
            "model_part_names" : [],
            "io_settings"  :{}
        }""")
        self.model = model
        self.parameters = parameters
        self.optimization_problem = optimization_problem
        self.parameters.ValidateAndAssignDefaults(default_settings)
        self.io_settings = self.parameters["io_settings"]
        self.s_connection_name = ""
        self.model_parts: 'list[Kratos.ModelPart]' = []

        for model_part_name in self.parameters["model_part_names"].GetStringArray():
            self.model.CreateModelPart(model_part_name)

    def GetAnalysisModelPart(self):
        return None

    @time_decorator()
    def Initialize(self):
        # Set up CoSim Connection Settings
        settings = CoSimIO.Info()
        settings.SetString("my_name", self.io_settings["my_name"].GetString())
        settings.SetString("connect_to", self.io_settings["connect_to"].GetString())
        settings.SetInt("echo_level", 1)
        settings.SetString("version", "1.25")
        settings.SetString("communication_format", self.io_settings["communication_format"].GetString())

        # Connecting
        return_info = CoSimIO.Connect(settings)
        self.s_connection_name = return_info.GetString("connection_name")

        # Reciving Design surface
        for model_part_name in self.parameters["model_part_names"].GetStringArray():
            # Send the signal
            self.Execute("ExportMesh")

            # Get the mesh
            info = CoSimIO.Info()
            info.SetString("connection_name", self.s_connection_name)
            info.SetString("identifier", model_part_name.replace(".", "-"))

            CoSimIO.ImportMesh(info, self.model[model_part_name], Kratos.ParallelEnvironment.GetDefaultDataCommunicator()) # doesn't support mpi, see cosim

        # initialize model parts
        self.model_parts = [self.model[model_part_name] for model_part_name in self.parameters["model_part_names"].GetStringArray()]

    def Check(self) -> None:
        pass

    def Finalize(self) -> None:
        self.analysis.Finalize()

    def Execute(self, command):
        # Send Current Shape
        # Solve Primal
        info = CoSimIO.Info()
        info.SetString("connection_name", self.s_connection_name)
        info.SetString("identifier", "run_control")
        info.SetString("control_signal", command)
        CoSimIO.ExportInfo(info)

    @time_decorator()
    def GetResponseValue(self, name):

        self.Execute("SolveSolutionStep")

        info = CoSimIO.Info()
        info.SetString("connection_name", self.s_connection_name)
        info.SetString("identifier", "info_exchange")
        imported_info = CoSimIO.ImportInfo(info)
        data = imported_info.GetDouble("value")
        print(f"{name}: {data}")
        return data

    def GetResponseGradients(self, name):

        self.Execute("SolveAdjointSolutionStep")

        for model_part in self.model_parts:
            info = CoSimIO.Info()
            info.SetString("connection_name", self.s_connection_name)
            info.SetString("identifier", "info_exchange")
            imported_data = CoSimIO.DoubleVector()
            CoSimIO.ImportData(info, imported_data)
            # print(imported_data)

    @staticmethod
    def __GetVariablesList(variable_names_list: 'list[str]') -> 'list[Any]':
        return [Kratos.KratosGlobals.GetVariable(variable_name) for variable_name in variable_names_list]



