# Importing the Kratos Library
import KratosMultiphysics as KM

# Importing the base class
from KratosMultiphysics.CoSimulationApplication.base_classes.co_simulation_solver_wrapper import CoSimulationSolverWrapper

# Other imports
from KratosMultiphysics.CoSimulationApplication.utilities import model_part_utilities
from KratosMultiphysics.CoSimulationApplication.utilities.data_communicator_utilities import GetRankZeroDataCommunicator

def Create(settings, model, solver_name):
    return FLOWerWrapper(settings, model, solver_name)

class FLOWerWrapper(CoSimulationSolverWrapper):
    """This class serves as wrapper for the CFD solver FLOWer
    """
    def __init__(self, settings, model, solver_name):
        super().__init__(settings, model, solver_name)

        settings_defaults = KM.Parameters("""{
            "model_parts_read"      : { },
            "model_parts_send"      : { },
            "model_parts_recv"      : { },
            "export_data"           : [ ],
            "import_data"           : [ ],
            "write_received_meshes" : false
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

            if self.settings["solver_wrapper_settings"]["write_received_meshes"].GetBool():
                KM.ModelPartIO(model_part_name, KM.IO.WRITE | KM.ModelPartIO.SKIP_TIMER).WriteModelPart(self.model[model_part_name])


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

    def AdvanceInTime(self, current_time):
        return 0.0 # TODO find a better solution here... maybe get time from solver through IO

    def _GetIOType(self):
        return "empire_io" # FLOWer currently only supports the EmpireIO

    def _GetDataCommunicator(self):
        # this solver does not support MPI
        # more specifically the EmpireIO does not support MPI
        # since FLOWer only uses the EmpireIO this has to be hardcoded (for now)
        return GetRankZeroDataCommunicator()
