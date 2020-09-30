# Importing the Kratos Library
import KratosMultiphysics as KM

# Importing the base class
from KratosMultiphysics.CoSimulationApplication.base_classes.co_simulation_solver_wrapper import CoSimulationSolverWrapper

# Other imports
from KratosMultiphysics.CoSimulationApplication.utilities import model_part_utilities

def Create(settings, model, solver_name):
    return FLOWerWrapper(settings, model, solver_name)

class FLOWerWrapper(CoSimulationSolverWrapper):
    """This class serves as wrapper for the CFD solver FLOWer
    """
    def __init__(self, settings, model, solver_name):
        super().__init__(settings, model, solver_name)

        settings_defaults = KM.Parameters("""{
            "model_parts_read" : { },
            "model_parts_send" : { },
            "model_parts_recv" : { }
        }""")

        self.settings["solver_wrapper_settings"].ValidateAndAssignDefaults(settings_defaults)

        model_part_utilities.CreateMainModelPartsFromCouplingDataSettings(self.settings["data"], self.model, self.name)
        model_part_utilities.AllocateHistoricalVariablesFromCouplingDataSettings(self.settings["data"], self.model, self.name)

    def Initialize(self):
        super().Initialize()

        for main_model_part_name, mdpa_file_name in self.settings["solver_wrapper_settings"]["model_parts_read"].items():
            KM.ModelPartIO(mdpa_file_name.GetString()).ReadModelPart(self.model[main_model_part_name])

        for model_part_name, comm_name in self.settings["solver_wrapper_settings"]["model_parts_send"].items():
            interface_config = {
                "comm_name" : comm_name.GetString(),
                "model_part_name" : model_part_name
            }
            self.ExportCouplingInterface(interface_config)

        for model_part_name, comm_name in self.settings["solver_wrapper_settings"]["model_parts_recv"].items():
            interface_config = {
                "comm_name" : comm_name.GetString(),
                "model_part_name" : model_part_name
            }

            self.ImportCouplingInterface(interface_config)

    def AdvanceInTime(self, current_time):
        return 0.0 # TODO find a better solution here... maybe get time from solver through IO

    def _GetIOType(self):
        return self.settings["io_settings"]["type"].GetString()