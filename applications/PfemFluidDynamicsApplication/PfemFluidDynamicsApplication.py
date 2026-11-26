
from KratosMultiphysics import _ImportApplication

# Applications required
from KratosDelaunayMeshingApplication import *

# Application dependent names and paths
from KratosPfemFluidDynamicsApplication import *
application = KratosPfemFluidDynamicsApplication()
application_name = "KratosPfemFluidDynamicsApplication"

_ImportApplication(application, application_name)
