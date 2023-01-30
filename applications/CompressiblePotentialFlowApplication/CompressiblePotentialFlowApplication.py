
# Application dependent names and paths
from KratosMultiphysics import _ImportApplication
import KratosMultiphysics.FluidDynamicsApplication
from KratosCompressiblePotentialFlowApplication import *
application = KratosCompressiblePotentialFlowApplication()
application_name = "KratosCompressiblePotentialFlowApplication"

_ImportApplication(application, application_name)
