import KratosMultiphysics
from KratosMultiphysics.kratos_utilities import CheckIfApplicationsAvailable
from KratosMultiphysics.kratos_utilities import IssueDeprecationWarning

from KratosMultiphysics.CoSimulationApplication.solver_wrappers.kratos.mpm_dirichlet_wrapper import MPMDirichletWrapper

# Importing MPM
if not CheckIfApplicationsAvailable("MPMApplication"):
    raise ImportError("The MPMApplication is not available!")
import KratosMultiphysics.MPMApplication as KPM
from KratosMultiphysics.MPMApplication.mpm_analysis import MPMAnalysis

def Create(settings, model, solver_name):
    IssueDeprecationWarning('CoSimulationApplication:','"ParticleMechanicsDirichletWrapper" is deprecated and replaced by "MPMDirichletWrapper"')
    return MPMDirichletWrapper(settings, model, solver_name)
