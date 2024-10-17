import KratosMultiphysics as KM
from KratosMultiphysics.analysis_stage import AnalysisStage

default_data_comm = KM.ParallelEnvironment.GetDefaultDataCommunicator()
if default_data_comm.IsDistributed():
    from KratosMultiphysics.CoSimulationApplication import MPIExtension

from KratosMultiphysics.CoSimulationApplication import CoSimIO

import sys

def FlushAll():
    sys.stdout.flush()
    KM.Logger.Flush()

def CreateAnalysisStageWithCoSimIO(BaseAnalysisStage):
    # concept from https://stackoverflow.com/a/1334242

    class AnalysisStageWithCoSimIO(BaseAnalysisStage):
        '''Adds coupling functionality with CoSimIO to an AnalysisStage'''

        def __init__(self, model, parameters):
            if not issubclass(BaseAnalysisStage, AnalysisStage):
                raise Exception(f'Given baseclass "{BaseAnalysisStage}" does not inherit from AnalysisStage!')

            super().__init__(model, parameters)

            # To avoid many prints
            severity = KM.Logger.Severity.WARNING
            if self.echo_level > 0:
                severity = KM.Logger.Severity.INFO
            KM.Logger.GetDefaultOutput().SetSeverity(severity)

            if self.echo_level > 0: FlushAll()

        def Initialize(self):
            super().Initialize()

            self.co_sim_settings = self.project_parameters["co_sim_settings"]
            self.is_strong_coupling = self.co_sim_settings["is_strong_coupling"].GetBool()

            connection_settings = CoSimIO.InfoFromParameters(self.project_parameters["co_sim_settings"]["io_settings"])

            if default_data_comm.IsDistributed():
                info = MPIExtension.CoSimIO.ConnectMPI(connection_settings, default_data_comm)
            else:
                info = CoSimIO.Connect(connection_settings)

            if self.echo_level > 0: FlushAll()

            self.connection_name = info.GetString("connection_name")
            if info.GetInt("connection_status") != CoSimIO.ConnectionStatus.Connected:
                raise Exception("Connecting failed!")

            self.communication_settings = self.co_sim_settings["communication_settings"]

            # Exporting meshes to CoSimulation
            for model_part_name in self.communication_settings["export_meshes"].GetStringArray():
                info = CoSimIO.Info()
                info.SetString("connection_name", self.connection_name)
                info.SetString("identifier", model_part_name.replace(".", "-"))

                CoSimIO.ExportMesh(info, self.model[model_part_name])

            if self.echo_level > 0: FlushAll()

        def RunSolutionLoop(self):
            """This function executes the solution loop of the AnalysisStage
            It can be overridden by derived classes
            """
            while self.KeepAdvancingSolutionLoop():
                self.time = self._GetSolver().AdvanceInTime(self.time)
                self.InitializeSolutionStep()
                self._GetSolver().Predict()

                if self.is_strong_coupling:
                    repeat_time_step = True
                    while repeat_time_step:
                        self.__InnerLoop()
                        info = CoSimIO.Info()
                        info.SetString("connection_name", self.connection_name)
                        info.SetString("identifier", "repeat_time_step_info")
                        repeat_time_step_info = CoSimIO.ImportInfo(info)
                        repeat_time_step = repeat_time_step_info.GetBool("repeat_time_step")
                else:
                    self.__InnerLoop()

                self.FinalizeSolutionStep()
                self.OutputSolutionStep()

        def Finalize(self):
            super().Finalize()

            disconnect_settings = CoSimIO.Info()
            disconnect_settings.SetString("connection_name", self.connection_name)
            info = CoSimIO.Disconnect(disconnect_settings)
            if info.GetInt("connection_status") != CoSimIO.ConnectionStatus.Disconnected:
                raise Exception("Disconnecting failed!")

        def _GetSimulationName(self):
            return super()._GetSimulationName() + " with CoSimIO"

        def __InnerLoop(self):
            if self.echo_level > 0: FlushAll()

            # Import fields
            for field_settings in self.communication_settings["import_fields"]:
                identifier = field_settings["identifier"].GetString()
                model_part_name = field_settings["model_part_name"].GetString()
                model_part = self.model[model_part_name]
                variable_name = field_settings["variable_name"].GetString()
                variable = KM.KratosGlobals.GetVariable(variable_name)

                info = CoSimIO.Info()
                info.SetString("connection_name", self.connection_name)
                info.SetString("identifier", identifier)
                CoSimIO.ImportData(info, model_part, variable, KM.Globals.DataLocation.NodeHistorical)

            if self.echo_level > 0: FlushAll()

            self._GetSolver().SolveSolutionStep()

            # Export fields
            for field_settings in self.communication_settings["export_fields"]:
                identifier = field_settings["identifier"].GetString()
                model_part_name = field_settings["model_part_name"].GetString()
                model_part = self.model[model_part_name]
                variable_name = field_settings["variable_name"].GetString()
                variable = KM.KratosGlobals.GetVariable(variable_name)

                info = CoSimIO.Info()
                info.SetString("connection_name", self.connection_name)
                info.SetString("identifier", identifier)
                CoSimIO.ExportData(info, model_part, variable, KM.Globals.DataLocation.NodeHistorical)

    return AnalysisStageWithCoSimIO
