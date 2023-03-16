# Importing the Kratos Library
import KratosMultiphysics as KM
from KratosMultiphysics.kratos_utilities import CheckIfApplicationsAvailable

# Importing the base class
from KratosMultiphysics.CoSimulationApplication.solver_wrappers.kratos import kratos_base_wrapper

# Importing Rom
if not CheckIfApplicationsAvailable("RomApplication"):
    raise ImportError("The RomApplication is not available!")
import KratosMultiphysics.RomApplication.rom_testing_utilities as rom_testing_utilities

def Create(settings, model, solver_name):
    return RomWrapper(settings, model, solver_name)

class RomWrapper(kratos_base_wrapper.KratosBaseWrapper):
    """This class is the interface to the RomApplication of Kratos"""

    def _CreateAnalysisStage(self):
        return rom_testing_utilities.SetUpSimulationInstance(self.model, self.project_parameters)

    def _GetDataCommunicator(self):
        # unfortunately the ConDiff solvers are using the global parallelism, it cannot be changed
        # to run with less cores or in serial with the current design!
        return super()._GetDataCommunicator()