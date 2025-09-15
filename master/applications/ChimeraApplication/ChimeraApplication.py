
# Application dependent names and paths
from KratosMultiphysics import _ImportApplication
import KratosMultiphysics.FluidDynamicsApplication
from KratosChimeraApplication import *
application = KratosChimeraApplication()
application_name = "KratosChimeraApplication"

_ImportApplication(application, application_name)
