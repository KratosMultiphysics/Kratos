#makes KratosMultiphysics backward compatible with python 2.6 and 2.7
from __future__ import print_function, absolute_import, division

from KratosMultiphysics import _ImportApplication

# Applications requiered
from KratosDelaunayMeshingApplication import *

# Application dependent names and paths
from KratosPfemFluidDynamicsApplication import *
application = KratosPfemFluidDynamicsApplication()
application_name = "KratosPfemFluidDynamicsApplication"

_ImportApplication(application, application_name)