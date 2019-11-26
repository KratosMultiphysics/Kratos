from __future__ import absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

from sys import argv

import KratosMultiphysics as KM

from KratosMultiphysics.StructuralMechanicsApplication.structural_mechanics_analysis import StructuralMechanicsAnalysis
import KratosMultiphysics.CoSimulationApplication as KratosCoSim
import KratosMultiphysics.CoSimulationApplication.co_simulation_tools as cs_tools

class StructuralMechanicsAnalysisWithCoSimIO(StructuralMechanicsAnalysis):
    '''Main script for structural mechanics with CoSimIO'''

    def __init__(self, model, parameters):
        super(StructuralMechanicsAnalysisWithCoSimIO,self).__init__(model, parameters)

        # To avoid many prints
        if (self.echo_level == 0):
            KM.Logger.GetDefaultOutput().SetSeverity(KM.Logger.Severity.WARNING)
        else:
            KM.Logger.GetDefaultOutput().SetSeverity(KM.Logger.Severity.INFO)

    def Initialize(self):
        super(StructuralMechanicsAnalysisWithCoSimIO, self).Initialize()
        self.co_sim_settings = self.project_parameters["co_sim_settings"]
        self.connection_name = self.co_sim_settings["connection_name"].GetString()
        self.is_strong_coupling = self.co_sim_settings["is_strong_coupling"].GetBool()

        KratosCoSim.CoSimIO.Connect(self.connection_name, cs_tools.ParametersToStringDict(self.project_parameters["co_sim_settings"]["io_settings"]))

        self.communication_settings = self.co_sim_settings["communication_settings"]

        # Exporting meshes to CoSimulation
        for model_part_name in self.communication_settings["export_meshes"].GetStringArray():
            KratosCoSim.CoSimIO.ExportMesh(self.connection_name, model_part_name, self.model[model_part_name])

    def __InnerLoop(self):
        # Import fields
        for field_settings in self.communication_settings["import_fields"]:
            identifier = field_settings["identifier"].GetString()
            model_part_name = field_settings["model_part_name"].GetString()
            model_part = self.model[model_part_name]
            variable_name = field_settings["variable_name"].GetString()
            variable = KM.KratosGlobals.GetVariable(variable_name)

            KratosCoSim.CoSimIO.ImportData(self.connection_name, identifier, model_part, variable, KratosCoSim.CoSimIO.DataLocation.NodeHistorical)

        self._GetSolver().SolveSolutionStep()

        # Export fields
        for field_settings in self.communication_settings["export_fields"]:
            identifier = field_settings["identifier"].GetString()
            model_part_name = field_settings["model_part_name"].GetString()
            model_part = self.model[model_part_name]
            variable_name = field_settings["variable_name"].GetString()
            variable = KM.KratosGlobals.GetVariable(variable_name)

            KratosCoSim.CoSimIO.ExportData(self.connection_name, identifier, model_part, variable, KratosCoSim.CoSimIO.DataLocation.NodeHistorical)


    def RunSolutionLoop(self):
        """This function executes the solution loop of the AnalysisStage
        It can be overridden by derived classes
        """
        while self.KeepAdvancingSolutionLoop():
            self.time = self._GetSolver().AdvanceInTime(self.time)
            self.InitializeSolutionStep()
            self._GetSolver().Predict()

            if self.is_strong_coupling:
                is_converged = False
                while not is_converged:
                    self.__InnerLoop()
                    is_converged = KratosCoSim.CoSimIO.IsConverged(self.connection_name)
            else:
                self.__InnerLoop()

            self.FinalizeSolutionStep()
            self.OutputSolutionStep()

    def Finalize(self):
        super(StructuralMechanicsAnalysisWithCoSimIO, self).Finalize()

        KratosCoSim.CoSimIO.Disconnect(self.connection_name)

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