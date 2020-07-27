from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
from KratosMultiphysics.kratos_utilities import CheckIfApplicationsAvailable

# Importing the base class
from KratosMultiphysics.CoSimulationApplication.solver_wrappers.kratos import kratos_base_wrapper

# Importing FluidDynamics
if not CheckIfApplicationsAvailable("FluidDynamicsApplication"):
    raise ImportError("The FluidDynamicsApplication is not available!")
from KratosMultiphysics.FluidDynamicsApplication.fluid_dynamics_analysis import FluidDynamicsAnalysis

def Create(settings, model, solver_name):
    return FluidDynamicsWrapper(settings, model, solver_name)

class FluidDynamicsWrapper(kratos_base_wrapper.KratosBaseWrapper):
    """This class is the interface to the FluidDynamicsApplication of Kratos"""

    def _CreateAnalysisStage(self):
        return FluidDynamicsAnalysis(self.model, self.project_parameters)
