# Importing the Kratos Library
from KratosMultiphysics.kratos_utilities import CheckIfApplicationsAvailable
import importlib

# Importing the base class
from KratosMultiphysics.CoSimulationApplication.solver_wrappers.kratos import kratos_base_wrapper

# Importing Rom
if not CheckIfApplicationsAvailable("RomApplication"):
    raise ImportError("The RomApplication is not available!")
import KratosMultiphysics.RomApplication.rom_analysis as RomAnalysis

def Create(settings, model, solver_name):
    return RomWrapper(settings, model, solver_name)

class RomWrapper(kratos_base_wrapper.KratosBaseWrapper):
    """This class is the interface to the RomApplication of Kratos"""

    def _CreateAnalysisStage(self):
        # Get the parent simulation class
        analysis_stage_module_name = self.project_parameters["analysis_stage"].GetString()
        analysis_stage_class_name = analysis_stage_module_name.split('.')[-1]
        analysis_stage_class_name = ''.join(x.title() for x in analysis_stage_class_name.split('_'))

        analysis_stage_module = importlib.import_module(analysis_stage_module_name)
        analysis_stage_class = getattr(analysis_stage_module, analysis_stage_class_name)
        instance_factory = RomAnalysis.CreateRomAnalysisInstance
        return instance_factory(analysis_stage_class, self.model, self.project_parameters)

    def _GetDataCommunicator(self):
        # unfortunately the ConDiff solvers are using the global parallelism, it cannot be changed
        # to run with less cores or in serial with the current design!
        return super()._GetDataCommunicator()