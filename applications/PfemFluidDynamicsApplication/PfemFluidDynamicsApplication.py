
from KratosMultiphysics import _ImportApplication

# Applications requiered
from KratosDelaunayMeshingApplication import *

# Application dependent names and paths
from KratosPfemFluidDynamicsApplication import *
application = KratosPfemFluidDynamicsApplication()
application_name = "KratosPfemFluidDynamicsApplication"

_ImportApplication(application, application_name)
