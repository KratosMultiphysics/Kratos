# Importing the Kratos Library
from KratosMultiphysics.kratos_utilities import CheckIfApplicationsAvailable

# Importing the base class
from KratosMultiphysics.CoSimulationApplication.solver_wrappers.kratos import kratos_base_wrapper

# Importing FluidDynamics
if not CheckIfApplicationsAvailable("ConvectionDiffusionApplication"):
    raise ImportError("The ConvectionDiffusionApplication is not available!")
from KratosMultiphysics.ConvectionDiffusionApplication.convection_diffusion_analysis import ConvectionDiffusionAnalysis

def Create(settings, model, solver_name):
    return ConvectionDiffusionWrapper(settings, model, solver_name)

class ConvectionDiffusionWrapper(kratos_base_wrapper.KratosBaseWrapper):
    """This class is the interface to the ConvectionDiffusionApplication of Kratos"""

    def _CreateAnalysisStage(self):
        return ConvectionDiffusionAnalysis(self.model, self.project_parameters)
