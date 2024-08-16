# Importing the Kratos Library
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
        model_part_utilities.CreateMainModelPartsFromCouplingDataSettings(self.settings["data"], self.model, self.name)
        model_part_utilities.AllocateHistoricalVariablesFromCouplingDataSettings(self.settings["data"], self.model, self.name)

    def Initialize(self):
        for model_part_name in self.settings["solver_wrapper_settings"]["import_meshes"].GetStringArray():
            interface_config = {"model_part_name" : model_part_name}
            self.ImportCouplingInterface(interface_config)

        super().Initialize()

        for data in self.data_dict.values():
            data.GetModelPart().GetRootModelPart().SetBufferSize(2)

    def AdvanceInTime(self, current_time):
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
    
    def ExportData(self, data_config):
        if data_config["type"] == "repeat_time_step":
            pass
        else:
            super().ExportData(data_config)
