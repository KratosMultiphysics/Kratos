# Importing the Kratos Library
import KratosMultiphysics as KM
from KratosMultiphysics.kratos_utilities import CheckIfApplicationsAvailable

# Importing the base class
from KratosMultiphysics.CoSimulationApplication.solver_wrappers.kratos import kratos_base_wrapper

# Using capwords()
import string

# The import module is used to import the AnalysisStage from the corresponding application
from importlib import import_module

def Create(settings, model, solver_name):
    # We assume that the location of analysis stage is in the same application as the solver or it is completely specified
    module_name = settings["analysis_stage"].GetString()
    imported_module = import_module(module_name)

    if settings.Has("analysis_name"):
        analysis_stage_name = settings["analysis_name"].GetString()
    else:
        # We assume that the name of the AnalysisStage is the same as the name of the module in PascalCase instead of cammel_case
        file_name = module_name.split(".")[-1]
        # Convert Snake case to Pascal case
        analysis_stage_name = string.capwords(file_name.replace("_", " ")).replace(" ", "")
    analysis = getattr(imported_module, analysis_stage_name)
    return GenericWrapper(settings, model, solver_name, analysis)

class GenericWrapper(kratos_base_wrapper.KratosBaseWrapper):
    """This class is the interface to a generic analysis of Kratos"""

    def __init__(self, settings, model, solver_name, analysis):
        self.analysis = analysis
        super().__init__(settings, model, solver_name)

    def _CreateAnalysisStage(self):
        return self.analysis(self.model, self.project_parameters)
