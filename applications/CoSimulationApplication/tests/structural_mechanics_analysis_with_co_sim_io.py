from sys import argv

import KratosMultiphysics as KM

from KratosMultiphysics.StructuralMechanicsApplication.structural_mechanics_analysis import StructuralMechanicsAnalysis
from KratosMultiphysics.CoSimulationApplication import CoSimIO

class StructuralMechanicsAnalysisWithCoSimIO(StructuralMechanicsAnalysis):
    '''Main script for structural mechanics with CoSimIO'''

    def __init__(self, model, parameters):
        super().__init__(model, parameters)

        # To avoid many prints
        if (self.echo_level == 0):
            KM.Logger.GetDefaultOutput().SetSeverity(KM.Logger.Severity.WARNING)
        else:
            KM.Logger.GetDefaultOutput().SetSeverity(KM.Logger.Severity.INFO)

    def Initialize(self):
        super().Initialize()
        self.co_sim_settings = self.project_parameters["co_sim_settings"]
        self.is_strong_coupling = self.co_sim_settings["is_strong_coupling"].GetBool()

        connection_settings = CoSimIO.InfoFromParameters(self.project_parameters["co_sim_settings"]["io_settings"])

        info = CoSimIO.Connect(connection_settings)
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

    def __InnerLoop(self):
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
            CoSimIO.ImportData(info, model_part, variable, CoSimIO.DataLocation.NodeHistorical)

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
            CoSimIO.ExportData(info, model_part, variable, CoSimIO.DataLocation.NodeHistorical)


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
        return "Structural Mechanics Analysis with CoSimIO"

if __name__ == '__main__':
    if len(argv) != 2:
        err_msg  = 'Wrong number of input arguments!\n'
        err_msg += 'Use this script in the following way:\n'
        err_msg += '    "python structural_mechanics_analysis_with_co_sim_io.py <project-parameter-file>.json"\n'
        raise Exception(err_msg)

    parameter_file_name = argv[1]

    with open(parameter_file_name,'r') as parameter_file:
        parameters = KM.Parameters(parameter_file.read())

    model = KM.Model()
    StructuralMechanicsAnalysisWithCoSimIO(model, parameters).Run()