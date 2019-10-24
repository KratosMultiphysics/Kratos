#makes KratosMultiphysics backward compatible with python 2.6 and 2.7
from __future__ import print_function, absolute_import, division

import KratosMultiphysics as KM

# Applications requiered
from KratosDelaunayMeshingApplication import *

# Application dependent names and paths
from KratosPfemFluidDynamicsApplication import *
application = KratosPfemFluidDynamicsApplication()
application_name = "KratosPfemFluidDynamicsApplication"
application_folder = "PfemFluidDynamicsApplication"

KM._ImportApplicationAsModule(application, application_name, application_folder, __path__)