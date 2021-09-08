# Importing the Kratos Library
from sys import argv
import KratosMultiphysics as KM

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

        self.settings["solver_wrapper_settings"].ValidateAndAssignDefaults(settings_defaults)
        model_part_utilities.CreateModelPartsFromCouplingDataSettings(self.settings["data"], self.model, self.name)
        model_part_utilities.AllocateHistoricalVariablesFromCouplingDataSettings(self.settings["data"], self.model, self.name)

    def Initialize(self):
        super(OpenFOAMWrapper, self).Initialize()

        for data in self.data_dict.values():
            data.GetModelPart().SetBufferSize(2)

        """ImportCouplingInterface()= Imports coupling interface from an external solver
        External solver sends, CoSimulation receives = load values"""
        for model_part_name in self.settings["solver_wrapper_settings"]["import_meshes"].GetStringArray():
            interface_config = {"model_part_name" : model_part_name}
            self.ImportCouplingInterface(interface_config)

    def AdvanceInTime(self, current_time):
        """What is the use of this method?"""
        return 0.0 # TODO find a better solution here... maybe get time from solver through IO

    def SolveSolutionStep(self):
        for data_name in self.settings["solver_wrapper_settings"]["import_data"].GetStringArray():
            data_config = {
                "type" : "coupling_interface_data",
                "interface_data" : self.GetInterfaceData(data_name)
            }
            self.ImportData(data_config)

        super().SolveSolutionStep()

        for data_name in self.settings["solver_wrapper_settings"]["export_data"].GetStringArray():
            data_config = {
                "type" : "coupling_interface_data",
                "interface_data" : self.GetInterfaceData(data_name)
            }
            self.ExportData(data_config)


    def _GetIOType(self):
        return self.settings["io_settings"]["type"].GetString()

    def Finalize(self):
        super().Finalize()

        """ # DisConnect from OpenFOAM adapter
        disconnect_settings = CoSimIO.Info()
        disconnect_settings.SetString("connection_name", self.connection_name)
        info = CoSimIO.Disconnect(disconnect_settings)
        if info.GetInt("connection_status") != CoSimIO.ConnectionStatus.Disconnected:
            raise Exception("Disconnecting failed!") """