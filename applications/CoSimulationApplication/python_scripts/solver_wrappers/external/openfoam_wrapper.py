# Importing the Kratos Library
import KratosMultiphysics as KM
import time
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
        print("\nKRATOS: __init__\n")
        super().__init__(settings, model, solver_name)

        settings_defaults = KM.Parameters("""{
            "import_meshes"            : [ ],
            "export_data"              : [ ],
            "import_data"              : [ ],
            "start_time"               : 0.0,
            "time_step"                : 0.0,
            "end_time"                 : 0.0,
            "strong_coupling"          : true,
            "first_solver_in_sequence" : false
                                          
        }""")

        self.settings["solver_wrapper_settings"].ValidateAndAssignDefaults(settings_defaults)
        model_part_utilities.CreateMainModelPartsFromCouplingDataSettings(self.settings["data"], self.model, self.name)
        model_part_utilities.AllocateHistoricalVariablesFromCouplingDataSettings(self.settings["data"], self.model, self.name)

        self.current_time = settings["solver_wrapper_settings"]["start_time"].GetDouble()
        self.time_step = settings["solver_wrapper_settings"]["time_step"].GetDouble()
        self.end_time = settings["solver_wrapper_settings"]["end_time"].GetDouble()
        self.isStrongCoupling = settings["solver_wrapper_settings"]["strong_coupling"].GetBool()
        self.isFirstSolverInSequence = settings["solver_wrapper_settings"]["first_solver_in_sequence"].GetBool()
        
        self.firstIteration = True

        print("\nKRATOS: __init__ Ende\n")

    def Initialize(self):
        print("\nKRATOS: Initialize\n")

        # Import meshes
        for model_part_name in self.settings["solver_wrapper_settings"]["import_meshes"].GetStringArray():
            interface_config = {"model_part_name" : model_part_name}
            self.ImportCouplingInterface(interface_config)

        # Test output modelPart
        # modelPart = self.model["interface"]
        # vtu_output: KM.VtuOutput = KM.VtkOutput(modelPart)
        # vtu_output.PrintOutput("test")

        super().Initialize()
        for data in self.data_dict.values():
            data.GetModelPart().GetRootModelPart().SetBufferSize(2)

        self.sendControlSignal("isStrongCoupling", {"isStrongCoupling" : self.isStrongCoupling})

        # Export time step settings
        self.sendControlSignal("setEndOfStepWindow", {"end_time_of_step_window": self.time_step})
        self.sendControlSignal("updateMaxTimeStepSize", {"maxTimeStepSize": self.time_step})
        self.sendControlSignal("firstOneToGo", {"firstOneToGo": self.isFirstSolverInSequence})

        print("\nKRATOS: Initialize Ende\n")

    def AdvanceInTime(self, current_time):
        print("AdvanceInTime!!!")

        if not self.firstIteration:
            self.current_time += self.time_step # adaptive time stepping might be different
        return self.current_time

    def sendControlSignal(self, command: str, additionalInfo: dict = {}):
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
                print("Datatype " + type(additionalInfo[key]) + " could not been sent.")
                print("Modify function \"sendControlSignal()\" in the adapter to add required datatype.")

        data_config = {
            "type": "control_signal",
            "control_signal": command,
            "settings": settings
        }
        self.ExportData(data_config)
        
    def SolveSolutionStep(self):
        print("\nKRATOS: SolveSolutionStep\n")
        print("  Solvername: " + str(self.settings["io_settings"]["co_sim_io_settings"]["my_name"].GetString()))

        # Import Data from first solver
        if self.firstIteration and self.isFirstSolverInSequence:
            print("    Importing Data " + self.settings["solver_wrapper_settings"]["import_data"].GetStringArray()[0])
            print("       from solver " + self.settings["io_settings"]["co_sim_io_settings"]["my_name"].GetString())

            for data_name in self.settings["solver_wrapper_settings"]["import_data"].GetStringArray():
                data_config = {
                "type" : "coupling_interface_data",
                "interface_data" : self.GetInterfaceData(data_name)
                }
                self.ImportData(data_config)

            self.firstIteration = False
            return
        
        elif self.firstIteration and not self.isFirstSolverInSequence:
            for data_name in self.settings["solver_wrapper_settings"]["export_data"].GetStringArray():
                data_config = {
                    "type" : "coupling_interface_data",
                    "interface_data" : self.GetInterfaceData(data_name)
                }
                print("    Export Data " + data_name)
                print("             To " + self.settings["io_settings"]["co_sim_io_settings"]["my_name"].GetString())
                self.ExportData(data_config)
                          
            for data_name in self.settings["solver_wrapper_settings"]["import_data"].GetStringArray():
                data_config = {
                "type" : "coupling_interface_data",
                "interface_data" : self.GetInterfaceData(data_name)
                }
                print("    Importing Data " + data_name)
                print("       from solver " + self.settings["io_settings"]["co_sim_io_settings"]["my_name"].GetString())
                self.ImportData(data_config)

            self.firstIteration = False
            return
        
        # Export the end time of the current time window
        self.sendControlSignal("setEndOfStepWindow", {"end_time_of_step_window": self.current_time + self.time_step})
        # print("Und updateMaxTimeStepSize()")
        # self.sendControlSignal("updateMaxTimeStepSize", {"maxTimeStepSize": self.time_step})


        # Export data to solver and solve
        print("Export data to solver and solve!")

        self.sendControlSignal("solve")
        for data_name in self.settings["solver_wrapper_settings"]["export_data"].GetStringArray():
            data_config = {
                "type" : "coupling_interface_data",
                "interface_data" : self.GetInterfaceData(data_name)
            }
            print("    Export Data " + data_name)
            print("             To " + self.settings["io_settings"]["co_sim_io_settings"]["my_name"].GetString())
            self.ExportData(data_config)

        print("Import Data")
        for data_name in self.settings["solver_wrapper_settings"]["import_data"].GetStringArray():
            data_config = {
                "type" : "coupling_interface_data",
                "interface_data" : self.GetInterfaceData(data_name)
            }
            print("    Importing Data " + data_name)
            print("       from solver " + self.settings["io_settings"]["co_sim_io_settings"]["my_name"].GetString())
            self.ImportData(data_config)
        
        print("\nKRATOS: SolveSolutionStep Ende\n")

    def Predict(self):
        print("Kratos: Predict")

    def InitializeSolutionStep(self):
        print("Kratos: InitializeSolutionStep")

    def FinalizeSolutionStep(self):
        print("Kratos: FinalizeSolutionStep")

    def OutputSolutionStep(self):
        print("Kratos: OutputSolutionStep")

    def Finalize(self):
        print("Kratos: Finalize")
        self.sendControlSignal("finalize")
        super().Finalize()


    def _GetIOType(self):
        return self.settings["io_settings"]["type"].GetString()
