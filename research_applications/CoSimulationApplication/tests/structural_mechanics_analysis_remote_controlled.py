from sys import argv

import KratosMultiphysics as KM

from KratosMultiphysics.StructuralMechanicsApplication.structural_mechanics_analysis import StructuralMechanicsAnalysis
from KratosMultiphysics.CoSimulationApplication import CoSimIO

from itertools import chain

class StructuralMechanicsAnalysisRemoteControlled(StructuralMechanicsAnalysis):
    '''Main script for structural mechanics remote controlled with CoSimIO'''

    def __init__(self, model, parameters):
        super().__init__(model, parameters)

        # To avoid many prints
        if self.echo_level == 0:
            KM.Logger.GetDefaultOutput().SetSeverity(KM.Logger.Severity.WARNING)
        else:
            KM.Logger.GetDefaultOutput().SetSeverity(KM.Logger.Severity.INFO)

    def Initialize(self):
        super().Initialize()

        connection_settings = CoSimIO.InfoFromParameters(self.project_parameters["co_sim_settings"]["io_settings"])

        info = CoSimIO.Connect(connection_settings)
        self.connection_name = info.GetString("connection_name")
        if info.GetInt("connection_status") != CoSimIO.ConnectionStatus.Connected:
            raise Exception("Connecting failed!")

        ## specifying which fields belong to which identifier when im-/exporting data
        self.data_comm_settings = {}
        for field_settings in chain(self.project_parameters["co_sim_settings"]["communication_settings"]["export_fields"],self.project_parameters["co_sim_settings"]["communication_settings"]["import_fields"]):
            identifier = field_settings["identifier"].GetString()
            model_part_name = field_settings["model_part_name"].GetString()
            model_part = self.model[model_part_name]
            variable_name = field_settings["variable_name"].GetString()
            variable = KM.KratosGlobals.GetVariable(variable_name)
            self.data_comm_settings[identifier] = {
                "model_part" : model_part,
                "variable" : variable
            }

        # declaring in place the functions that are registered by the CoSimIO
        def _ImportData(info):
            identifier = info.GetString("identifier")
            model_part = self.data_comm_settings[identifier]["model_part"]
            variable   = self.data_comm_settings[identifier]["variable"]

            info = CoSimIO.Info()
            info.SetString("connection_name", self.connection_name)
            info.SetString("identifier", identifier)
            CoSimIO.ImportData(info, model_part, variable, CoSimIO.DataLocation.NodeHistorical)
            return CoSimIO.Info()

        def _ExportData(info):
            identifier = info.GetString("identifier")
            model_part = self.data_comm_settings[identifier]["model_part"]
            variable   = self.data_comm_settings[identifier]["variable"]

            info = CoSimIO.Info()
            info.SetString("connection_name", self.connection_name)
            info.SetString("identifier", identifier)
            CoSimIO.ExportData(info, model_part, variable, CoSimIO.DataLocation.NodeHistorical)
            return CoSimIO.Info()

        def _ExportMesh(info):
            model_part_name = info.GetString("identifier")

            export_info = CoSimIO.Info()
            export_info.SetString("connection_name", self.connection_name)
            export_info.SetString("identifier", model_part_name.replace(".", "-"))

            CoSimIO.ExportMesh(export_info, self.model[model_part_name])
            return CoSimIO.Info()

        def _AdvanceInTime(info):
            current_time = info.GetDouble("current_time")
            self.time = self._GetSolver().AdvanceInTime(current_time)
            return CoSimIO.Info()

        def _InitializeSolutionStep(info):
            self.InitializeSolutionStep()
            return CoSimIO.Info()

        def _Predict(info):
            self._GetSolver().Predict()
            return CoSimIO.Info()

        def _SolveSolutionStep(info):
            self._GetSolver().SolveSolutionStep()
            return CoSimIO.Info()

        def _FinalizeSolutionStep(info):
            self.FinalizeSolutionStep()
            return CoSimIO.Info()

        def _OutputSolutionStep(info):
            self.OutputSolutionStep()
            return CoSimIO.Info()

        register_info = CoSimIO.Info()
        register_info.SetString("connection_name", self.connection_name)

        register_info.SetString("function_name", "ImportData")
        CoSimIO.Register(register_info, _ImportData)

        register_info.SetString("function_name", "ExportData")
        CoSimIO.Register(register_info, _ExportData)

        register_info.SetString("function_name", "ExportMesh")
        CoSimIO.Register(register_info, _ExportMesh)

        register_info.SetString("function_name", "AdvanceInTime")
        CoSimIO.Register(register_info, _AdvanceInTime)

        register_info.SetString("function_name", "InitializeSolutionStep")
        CoSimIO.Register(register_info, _InitializeSolutionStep)

        register_info.SetString("function_name", "Predict")
        CoSimIO.Register(register_info, _Predict)

        register_info.SetString("function_name", "SolveSolutionStep")
        CoSimIO.Register(register_info, _SolveSolutionStep)

        register_info.SetString("function_name", "FinalizeSolutionStep")
        CoSimIO.Register(register_info, _FinalizeSolutionStep)

        register_info.SetString("function_name", "OutputSolutionStep")
        CoSimIO.Register(register_info, _OutputSolutionStep)

    def RunSolutionLoop(self):
        # running remote controlled simulation
        # deliberately not calling baseclass as this is the remote controlled case
        run_info = CoSimIO.Info()
        run_info.SetString("connection_name", self.connection_name)
        CoSimIO.Run(run_info)

    def Finalize(self):
        super().Finalize()

        disconnect_settings = CoSimIO.Info()
        disconnect_settings.SetString("connection_name", self.connection_name)
        info = CoSimIO.Disconnect(disconnect_settings)
        if info.GetInt("connection_status") != CoSimIO.ConnectionStatus.Disconnected:
            raise Exception("Disconnecting failed!")

    def _GetSimulationName(self):
        return "Structural Mechanics Analysis remote controlled with CoSimIO"


if __name__ == '__main__':
    if len(argv) != 2:
        err_msg  = 'Wrong number of input arguments!\n'
        err_msg += 'Use this script in the following way:\n'
        err_msg += '    "python structural_mechanics_analysis_remote_controlled.py <project-parameter-file>.json"\n'
        raise Exception(err_msg)

    parameter_file_name = argv[1]

    with open(parameter_file_name,'r') as parameter_file:
        parameters = KM.Parameters(parameter_file.read())

    model = KM.Model()
    StructuralMechanicsAnalysisRemoteControlled(model, parameters).Run()
