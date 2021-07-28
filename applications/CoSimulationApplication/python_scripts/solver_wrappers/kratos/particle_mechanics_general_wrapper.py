# Importing the Kratos Library
from KratosMultiphysics.kratos_utilities import CheckIfApplicationsAvailable

# Importing the base class
from KratosMultiphysics.CoSimulationApplication.solver_wrappers.kratos import kratos_base_wrapper

# Importing ParticleMechanics
if not CheckIfApplicationsAvailable("ParticleMechanicsApplication"):
    raise ImportError("The ParticleMechanicsApplication is not available!")
from KratosMultiphysics.ParticleMechanicsApplication.particle_mechanics_analysis import ParticleMechanicsAnalysis

def Create(settings, model, solver_name):
    return ParticleMechanicsWrapper(settings, model, solver_name)

class ParticleMechanicsWrapper(kratos_base_wrapper.KratosBaseWrapper):
    """This class is the interface to the ParticleMechanicsApplication of Kratos"""

    def _CreateAnalysisStage(self):
        return ParticleMechanicsAnalysis(self.model, self.project_parameters)
