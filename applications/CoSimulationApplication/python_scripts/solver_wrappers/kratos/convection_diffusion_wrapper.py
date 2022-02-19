# Importing the Kratos Library
import KratosMultiphysics as KM
from KratosMultiphysics.kratos_utilities import CheckIfApplicationsAvailable

# Importing the base class
from KratosMultiphysics.CoSimulationApplication.solver_wrappers.kratos import kratos_base_wrapper

# Importing ConvectionDiffusion
if not CheckIfApplicationsAvailable("ConvectionDiffusionApplication"):
    raise ImportError("The ConvectionDiffusionApplication is not available!")
from KratosMultiphysics.ConvectionDiffusionApplication.convection_diffusion_analysis import ConvectionDiffusionAnalysis

def Create(settings, model, solver_name):
    return ConvectionDiffusionWrapper(settings, model, solver_name)

class ConvectionDiffusionWrapper(kratos_base_wrapper.KratosBaseWrapper):
    """This class is the interface to the ConvectionDiffusionApplication of Kratos"""

    def _CreateAnalysisStage(self):
        return ConvectionDiffusionAnalysis(self.model, self.project_parameters)

    def _GetDataCommunicator(self):
        # unfortunately the ConDiff solvers are using the global parallelism, it cannot be changed
        # to run with less cores or in serial with the current design!
        return super()._GetDataCommunicator()
