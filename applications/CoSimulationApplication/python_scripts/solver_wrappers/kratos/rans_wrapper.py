# Importing the Kratos Library
from KratosMultiphysics.kratos_utilities import CheckIfApplicationsAvailable

# Importing the base class
from KratosMultiphysics.CoSimulationApplication.solver_wrappers.kratos import kratos_base_wrapper

# Importing FluidDynamics
if not CheckIfApplicationsAvailable("RANSApplication"):
    raise ImportError("The RANSApplication is not available!")
from KratosMultiphysics.RANSApplication.rans_analysis import RANSAnalysis

def Create(settings, model, solver_name):
    return RANSWrapper(settings, model, solver_name)

class RANSWrapper(kratos_base_wrapper.KratosBaseWrapper):
    """This class is the interface to the RANSApplication of Kratos"""

    def _CreateAnalysisStage(self):
        return RANSAnalysis(self.model, self.project_parameters)
