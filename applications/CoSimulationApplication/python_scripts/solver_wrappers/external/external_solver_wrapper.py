# Importing the Kratos Library
import KratosMultiphysics as KM
import KratosMultiphysics.StructuralMechanicsApplication # needed for some variables

# Importing the base class
from KratosMultiphysics.CoSimulationApplication.base_classes.co_simulation_solver_wrapper import CoSimulationSolverWrapper

# Other imports
from KratosMultiphysics.CoSimulationApplication.utilities import model_part_utilities

def Create(settings, model, solver_name):
    return ExternalSolverWrapper(settings, model, solver_name)

class ExternalSolverWrapper(CoSimulationSolverWrapper):
    """This class is a generic wrapper for connecting external solvers
    The import of meshes is done once in the beginning
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