# Importing the Kratos Library
import KratosMultiphysics as Kratos
from KratosMultiphysics.CoSimulationApplication.coupling_interface_data import CouplingInterfaceData
import KratosMultiphysics.StructuralMechanicsApplication # needed for some variables
import importlib

# Importing the base class
from KratosMultiphysics.CoSimulationApplication.base_classes.co_simulation_solver_wrapper import CoSimulationSolverWrapper

# Other imports
from KratosMultiphysics.CoSimulationApplication.utilities import model_part_utilities

def Create(settings: Kratos.Parameters, model: Kratos.Model, model_part_name: str):
    return ExternalPythonSolverWrapper(settings, model, model_part_name)

class ExternalPythonSolverWrapper(CoSimulationSolverWrapper):
    """This class is a generic wrapper for connecting external solvers
    The import of meshes is done once in the beginning
    """
    def __init__(self, settings: Kratos.Parameters, model: Kratos.Model, model_part_name: str):
        super().__init__(settings, model, model_part_name)

        settings_defaults = Kratos.Parameters("""{
            "python_module"    : "Path_to_package",
            "model_part_name"  : "Single_node",
            "import_meshes"    : [ ],
            "export_data"      : [ ],
            "import_data"      : [ ]
        }""")


        self.settings["solver_wrapper_settings"].ValidateAndAssignDefaults(settings_defaults)

        package_path = self.settings["solver_wrapper_settings"]["python_module"].GetString()
        self.module = importlib.import_module(package_path)
        model_part_name = self.settings["solver_wrapper_settings"]["model_part_name"].GetString()
        self.model_part = self.model.CreateModelPart(model_part_name)
        model_part_utilities.CreateMainModelPartsFromCouplingDataSettings(self.settings["data"], self.model, self.name)
        model_part_utilities.AllocateHistoricalVariablesFromCouplingDataSettings(self.settings["data"], self.model, self.name)
        
        self.model_part.CreateNewNode(1, 0.0, 0.0, 0.0)

        self.external_model = self.module.CreateModel()

    def Initialize(self):
        super().Initialize()
        self.external_model.Initialize()

    def AdvanceInTime(self, current_time):
        return 1.0 # TODO find a better solution here... maybe get time from solver through IO

    def SolveSolutionStep(self):
        print("Start SolveSolutionStep")
        for data_name in self.settings["solver_wrapper_settings"]["export_data"].GetStringArray():
            data_config = {
                "type" : "direct_value",
                "external_model": self.external_model,
                "interface_data" : self.GetInterfaceData(data_name)
            }
            self.ExportData(data_config)

        self.external_model.SolveSolutionStep()
        super().SolveSolutionStep()

        for data_name in self.settings["solver_wrapper_settings"]["import_data"].GetStringArray():
            data_config = {
                "type" : "direct_value",
                "external_model": self.external_model,
                "interface_data" : self.GetInterfaceData(data_name)
            }
            self.ImportData(data_config)

    def _GetIOType(self):
        return self.settings["io_settings"]["type"].GetString()