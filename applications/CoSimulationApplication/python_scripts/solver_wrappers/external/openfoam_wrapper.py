# Importing the Kratos Library
from sys import argv
import KratosMultiphysics as KM
from KratosMultiphysics.CoSimulationApplication import CoSimIO

# Importing the base class
from KratosMultiphysics.CoSimulationApplication.base_classes.co_simulation_solver_wrapper import CoSimulationSolverWrapper

# Other imports
from KratosMultiphysics.CoSimulationApplication.utilities import model_part_utilities

def Create(settings, model, solver_name):
    return OpenFOAMWrapper(settings, model, solver_name)

class OpenFOAMWrapper(CoSimulationSolverWrapper):
    """This class serves as wrapper for the CFD solver OpenFOAM
    """
    def __init__(self, settings, model, solver_name):
        super().__init__(settings, model, solver_name)

        settings_defaults = KM.Parameters("""{
            "import_meshes"    : [ ],
            "export_data"      : [ ],
            "import_data"      : [ ]
        }""")

        #self.settings["solver_wrapper_settings"].ValidateAndAssignDefaults(settings_defaults)
        model_part_utilities.CreateModelPartsFromCouplingDataSettings(self.settings["data"], self.model, self.name)
        model_part_utilities.AllocateHistoricalVariablesFromCouplingDataSettings(self.settings["data"], self.model, self.name)

    def Initialize(self):
        super().Initialize()

        # Check which IO is used
        self.type_io_used =  self.project_parameters["openfoam_kratos_wrapper"]["io_settings"]["type"]

        # Connect to OpenFOAM adapter
        connection_settings = CoSimIO.InfoFromParameters(self.project_parameters["openfoam_kratos_wrapper"]["io_settings"])
        info = CoSimIO.Connect(connection_settings)
        self.connection_name = info.GetString("connection_name")
        if info.GetInt("connection_status") != CoSimIO.ConnectionStatus.Connected:
            raise Exception("Connecting failed!")

        for model_part_name in self.settings["solver_wrapper_settings"]["import_meshes"].GetStringArray():
            interface_config = { "model_part_name" : model_part_name }
            self.ImportCouplingInterface(interface_config)

        self.communication_settings = self.co_sim_settings["communication_settings"]

        # Importing meshes from OpenFOAM to CoSimulation
        for model_part_name in self.communication_settings["export_meshes"].GetStringArray():
            info = CoSimIO.Info()
            info.SetString("connection_name", self.connection_name)
            info.SetString("identifier", model_part_name.replace(".", "-"))

            CoSimIO.ImportData(info, self.model[model_part_name])

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


    def AdvanceInTime(self, current_time):
        return 0.0 # TODO find a better solution here... maybe get time from solver through IO

    def SolveSolutionStep(self):
        for data_name in self.settings["solver_wrapper_settings"]["export_data"].GetStringArray():
            data_config = {
                "type" : "coupling_interface_data",
                "interface_data" : self.GetInterfaceData(data_name)
            }
            self.ExportData(data_config)

        super().SolveSolutionStep()

        for data_name in self.settings["solver_wrapper_settings"]["import_data"].GetStringArray():
            data_config = {
                "type" : "coupling_interface_data",
                "interface_data" : self.GetInterfaceData(data_name)
            }
            self.ImportData(data_config)

    def _GetIOType(self):
        return self.settings["io_settings"]["type"].GetString()

    def Finalize(self):
        super().Finalize()

        # DisConnect from OpenFOAM adapter
        disconnect_settings = CoSimIO.Info()
        disconnect_settings.SetString("connection_name", self.connection_name)
        info = CoSimIO.Disconnect(disconnect_settings)
        if info.GetInt("connection_status") != CoSimIO.ConnectionStatus.Disconnected:
            raise Exception("Disconnecting failed!")
