# Importing the Kratos Library
from KratosMultiphysics.kratos_utilities import CheckIfApplicationsAvailable

# Importing the base class
from KratosMultiphysics.CoSimulationApplication.solver_wrappers.kratos import kratos_base_wrapper

# Importing PfemFluidDynamics
if not CheckIfApplicationsAvailable("PfemFluidDynamicsApplication"):
    raise ImportError("The PfemFluidDynamicsApplication is not available!")
from KratosMultiphysics.PfemFluidDynamicsApplication.pfem_fluid_dynamics_analysis import PfemFluidDynamicsAnalysis

def Create(settings, model, solver_name):
    return PfemFluidDynamicsWrapper(settings, model, solver_name)

class PfemFluidDynamicsWrapper(kratos_base_wrapper.KratosBaseWrapper):
    """This class is the interface to the PfemFluidDynamicsApplication of Kratos"""

    def _CreateAnalysisStage(self):
        return PfemFluidDynamicsAnalysis(self.model, self.project_parameters)