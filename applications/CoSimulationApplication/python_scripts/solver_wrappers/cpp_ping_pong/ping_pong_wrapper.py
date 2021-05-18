import subprocess
# Importing the Kratos Library
import KratosMultiphysics as KM

# Importing the base class
from KratosMultiphysics.CoSimulationApplication.base_classes.co_simulation_solver_wrapper import CoSimulationSolverWrapper

# Other imports
from KratosMultiphysics.CoSimulationApplication.utilities import model_part_utilities

def Create(settings, model, solver_name):
    return PingPongWrapper(settings, model, solver_name)

class PingPongWrapper(CoSimulationSolverWrapper):
    """This class serves as wrapper for the cpp ping and pong solvers
    """
    def __init__(self, settings, model, solver_name):
        super().__init__(settings, model, solver_name)

        settings_defaults = KM.Parameters("""{
            "main_model_part_name" : "",
            "domain_size" : 2,
            "executable_name"  : "",
            "export_data"      : [ ],
            "import_data"      : [ ]
        }""")

        self.settings["solver_wrapper_settings"].ValidateAndAssignDefaults(settings_defaults)
        model_part_name = self.settings["solver_wrapper_settings"]["main_model_part_name"].GetString()
        model_part_utilities.CreateMainModelPartsFromCouplingDataSettings(self.settings["data"], self.model, self.name)
        model_part_utilities.AllocateHistoricalVariablesFromCouplingDataSettings(self.settings["data"], self.model, self.name)
        self.mp = self.model[model_part_name]
        self.mp.CreateNewNode(1,0,0,0)

    def Initialize(self):
        super().Initialize()

    def Finalize(self):
        super().Finalize()
        with self.rv.stdout, open(self.name +'.log', 'w') as file:
            for line in self.rv.stdout:
                file.write(line.decode("utf-8"))

    def AdvanceInTime(self, current_time):
        return 1.0

    def SolveSolutionStep(self):
        for data_name in self.settings["solver_wrapper_settings"]["export_data"].GetStringArray():
            data_config = {
                "type" : "coupling_interface_data",
                "interface_data" : self.GetInterfaceData(data_name)
            }
            self.ExportData(data_config)

        super().SolveSolutionStep()
        self.__RunExecutable()

        for data_name in self.settings["solver_wrapper_settings"]["import_data"].GetStringArray():
            data_config = {
                "type" : "coupling_interface_data",
                "interface_data" : self.GetInterfaceData(data_name)
            }
            self.ImportData(data_config)

    def _GetIOType(self):
        return self.settings["io_settings"]["type"].GetString()

    def __RunExecutable(self):
        command_txt = self.settings["solver_wrapper_settings"]["executable_name"].GetString()
        self.rv = subprocess.Popen(command_txt, stdout = subprocess.PIPE, stderr = subprocess.PIPE, shell=True, start_new_session=True)
