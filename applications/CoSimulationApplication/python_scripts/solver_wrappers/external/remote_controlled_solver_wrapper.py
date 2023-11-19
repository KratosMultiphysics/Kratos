# CoSimulation imports
import KratosMultiphysics as KM
import KratosMultiphysics.CoSimulationApplication as KratosCoSim
import KratosMultiphysics.StructuralMechanicsApplication # needed for some variables

# Importing the base class
from KratosMultiphysics.CoSimulationApplication.base_classes.co_simulation_solver_wrapper import CoSimulationSolverWrapper

# Other imports
from KratosMultiphysics.CoSimulationApplication.utilities import model_part_utilities

def Create(settings, model, solver_name):
    return RemoteControlledSolverWrapper(settings, model, solver_name)

class RemoteControlledSolverWrapper(CoSimulationSolverWrapper):
    """This class is a generic wrapper for connecting external solvers that are being remote controlled
    """
    def __init__(self, settings, model, solver_name):
        super().__init__(settings, model, solver_name)

        settings_defaults = KM.Parameters("""{
            "import_meshes"    : [ ],
            "export_data"      : [ ],
            "import_data"      : [ ]
        }""")

        self.settings["solver_wrapper_settings"].ValidateAndAssignDefaults(settings_defaults)

        model_part_utilities.CreateModelPartsFromCouplingDataSettings(self.settings["data"], self.model, self.name)
        model_part_utilities.AllocateHistoricalVariablesFromCouplingDataSettings(self.settings["data"], self.model, self.name)

    def Initialize(self):
        super().Initialize()

        for model_part_name in self.settings["solver_wrapper_settings"]["import_meshes"].GetStringArray():
            interface_config = { "model_part_name" : model_part_name }
            self.ImportCouplingInterface(interface_config)

    def AdvanceInTime(self, current_time):
        settings = KM.Parameters("""{}""")
        settings.AddEmptyValue("current_time").SetDouble(current_time)
        self.__SendControlSignal("AdvanceInTime", settings)

        # here one could import back the new time from the solver
        return 0.0

    def Predict(self):
        super().Predict()
        self.__SendControlSignal("Predict")

    def InitializeSolutionStep(self):
        super().InitializeSolutionStep()
        self.__SendControlSignal("InitializeSolutionStep")

    def SolveSolutionStep(self):
        for data_name in self.settings["solver_wrapper_settings"]["export_data"].GetStringArray():
            data_config = {
                "type" : "coupling_interface_data",
                "interface_data" : self.GetInterfaceData(data_name)
            }
            self.ExportData(data_config)

        super().SolveSolutionStep()
        self.__SendControlSignal("SolveSolutionStep")

        for data_name in self.settings["solver_wrapper_settings"]["import_data"].GetStringArray():
            data_config = {
                "type" : "coupling_interface_data",
                "interface_data" : self.GetInterfaceData(data_name)
            }
            self.ImportData(data_config)

    def FinalizeSolutionStep(self):
        super().FinalizeSolutionStep()
        self.__SendControlSignal("FinalizeSolutionStep")

    def OutputSolutionStep(self):
        super().OutputSolutionStep()
        self.__SendControlSignal("OutputSolutionStep")

    def Finalize(self):
        self.__SendControlSignal("exit")
        super().Finalize() # this also does the disconnect

    def ImportCouplingInterface(self, interface_config):
        settings = KM.Parameters("""{}""")
        settings.AddEmptyValue("identifier").SetString(interface_config["model_part_name"])
        self.__SendControlSignal("ExportMesh", settings) # TODO this can also be geometry at some point
        super().ImportCouplingInterface(interface_config)

    def ExportCouplingInterface(self, interface_config):
        self.__SendControlSignal("ImportMesh", interface_config["model_part_name"]) # TODO this can also be geometry at some point
        super().ExportCouplingInterface(interface_config)

    def ImportData(self, data_config):
        # CoSim imports, the external solver exports
        settings = KM.Parameters("""{}""")
        settings.AddEmptyValue("identifier").SetString(data_config["interface_data"].name)
        self.__SendControlSignal("ExportData", settings)
        super().ImportData(data_config)

    def ExportData(self, data_config):
        if data_config["type"] == "coupling_interface_data":
            # CoSim exports, the external solver imports
            settings = KM.Parameters("""{}""")
            settings.AddEmptyValue("identifier").SetString(data_config["interface_data"].name)
            self.__SendControlSignal("ImportData", settings)
        elif data_config["type"] == "repeat_time_step":
            return # we control the ext solver, no need for sending a repeat_time_step signal
        super().ExportData(data_config)

    def _GetIOType(self):
        return "kratos_co_sim_io"

    def __SendControlSignal(self, signal, settings=None):
        data_config = {
            "type"           : "control_signal",
            "control_signal" : signal,
            "settings"       : settings
        }
        self.ExportData(data_config)
