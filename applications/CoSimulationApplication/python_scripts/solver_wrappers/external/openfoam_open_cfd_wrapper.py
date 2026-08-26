# Importing the Kratos Library
import KratosMultiphysics as KM
import time
# Importing the base class
from KratosMultiphysics.CoSimulationApplication.base_classes.co_simulation_solver_wrapper import CoSimulationSolverWrapper

# Other imports
from KratosMultiphysics.CoSimulationApplication.utilities import model_part_utilities

def Create(settings, model, solver_name):
    return OpenFOAMOpenCFDWrapper(settings, model, solver_name)

class OpenFOAMOpenCFDWrapper(CoSimulationSolverWrapper):
    """@brief Co-simulation wrapper for the OpenCFD distribution of OpenFOAM.

    This wrapper supports both weak and strong coupling and is intended to be
    used with the OpenFOAM CoSimIO adapter:
    https://github.com/juancamarotti/OpenFOAM_CoSimIO-Adapter

    @note The adapter was implemented using OpenFOAM v2512. Compatibility with
    other OpenFOAM versions is not guaranteed.
    """
    def __init__(self, settings, model, solver_name):
        super().__init__(settings, model, solver_name)

        settings_defaults = KM.Parameters("""{
            "import_meshes"            : [ ],
            "export_data"              : [ ],
            "import_data"              : [ ],
            "start_time"               : 0.0,
            "time_step"                : 0.0,
            "end_time"                 : 0.0,
            "strong_coupling"          : true               
        }""")

        solver_wrapper_settings = self.settings["solver_wrapper_settings"]
        if not solver_wrapper_settings.Has("time_step"):
            raise ValueError(
                '"time_step" must be specified in "solver_wrapper_settings" '
                'for the OpenFOAM OpenCFD wrapper.'
            )

        solver_wrapper_settings.ValidateAndAssignDefaults(settings_defaults)
        model_part_utilities.CreateMainModelPartsFromCouplingDataSettings(self.settings["data"], self.model, self.name)
        model_part_utilities.AllocateHistoricalVariablesFromCouplingDataSettings(self.settings["data"], self.model, self.name)

        self.current_time = solver_wrapper_settings["start_time"].GetDouble()
        self.time_step = solver_wrapper_settings["time_step"].GetDouble()
        self.end_time = solver_wrapper_settings["end_time"].GetDouble()
        self.is_strong_coupling = solver_wrapper_settings["strong_coupling"].GetBool()

        if self.time_step <= 0.0:
            raise ValueError(
                '"time_step" must be greater than zero for the OpenFOAM '
                'OpenCFD wrapper.'
            )
        
        self.first_iteration = True

    def Initialize(self):
        # Import meshes
        for model_part_name in self.settings["solver_wrapper_settings"]["import_meshes"].GetStringArray():
            interface_config = {"model_part_name" : model_part_name}
            self.ImportCouplingInterface(interface_config)

        super().Initialize()
        for data in self.data_dict.values():
            data.GetModelPart().GetRootModelPart().SetBufferSize(2)

        # Import initial coupling data from OpenFOAM
        for data_name in self.settings["solver_wrapper_settings"]["import_data"].GetStringArray():
            data_config = {
                "type" : "coupling_interface_data",
                "interface_data" : self.GetInterfaceData(data_name)
            }
            self.ImportData(data_config)

        # Send a signal to define if the coupling is strong or not
        self.SendControlSignal("is_strong_coupling", {
            "is_strong_coupling": self.is_strong_coupling
        })

        # Export initial coupling data to OpenFOAM
        for data_name in self.settings["solver_wrapper_settings"]["export_data"].GetStringArray():
            data_config = {
                "type" : "coupling_interface_data",
                "interface_data" : self.GetInterfaceData(data_name)
            }
            self.ExportData(data_config)

    def AdvanceInTime(self, current_time):
        self.current_time = current_time + self.time_step
        return self.current_time

    def SendControlSignal(self, command: str, additionalInfo: dict = {}):
        keys = list(additionalInfo.keys())

        settings = KM.Parameters("""{}""")
        for key in keys:
            if type(additionalInfo[key]) == str:
                settings.AddEmptyValue(key).SetString(additionalInfo[key])
            elif type(additionalInfo[key]) == bool:
                settings.AddEmptyValue(key).SetBool(additionalInfo[key])
            elif type(additionalInfo[key]) == float:
                settings.AddEmptyValue(key).SetDouble(additionalInfo[key])
            elif type(additionalInfo[key]) == int:
                settings.AddEmptyValue(key).SetInt(additionalInfo[key])
            else:
                pass

        data_config = {
            "type": "control_signal",
            "control_signal": command,
            "settings": settings
        }
        self.ExportData(data_config)
        
    def SolveSolutionStep(self):
        # In the first coupling iteration we only import data
        if self.first_iteration:
            for data_name in self.settings["solver_wrapper_settings"]["import_data"].GetStringArray():
                data_config = {
                "type" : "coupling_interface_data",
                "interface_data" : self.GetInterfaceData(data_name)
                }
                self.ImportData(data_config)

            self.first_iteration = False
            return

        # Export data to the adapter
        for data_name in self.settings["solver_wrapper_settings"]["export_data"].GetStringArray():
            data_config = {
                "type" : "coupling_interface_data",
                "interface_data" : self.GetInterfaceData(data_name)
            }
            self.ExportData(data_config)

        # Import data from the adapter
        for data_name in self.settings["solver_wrapper_settings"]["import_data"].GetStringArray():
            data_config = {
                "type" : "coupling_interface_data",
                "interface_data" : self.GetInterfaceData(data_name)
            }
            self.ImportData(data_config)

    def FinalizeSolutionStep(self):
        tol = 1e-10

        is_last_step = self.current_time + tol >= self.end_time

        self.SendControlSignal(
            "coupling_ongoing",
            {"coupling_ongoing": not is_last_step}
        )

        super().FinalizeSolutionStep()

    def _GetIOType(self):
        return self.settings["io_settings"]["type"].GetString()
    
    def ExportData(self, data_config):
        if data_config["type"] == "repeat_time_step":
            self._SendRepeatTimeStep(data_config["repeat_time_step"])
            return

        super().ExportData(data_config)

    # This function communicates to the external solver whether the time step needs to be repeated or not
    def _SendRepeatTimeStep(self, repeat_time_step):
        settings = KM.Parameters("""{}""")
        settings.AddEmptyValue("repeat_time_step").SetBool(repeat_time_step)

        self.SendControlSignal(
            "repeat_time_step",
            {"repeat_time_step": repeat_time_step}
        )
