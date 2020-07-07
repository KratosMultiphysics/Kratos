# makes KratosMultiphysics backward compatible with python 2.6 and 2.7
from __future__ import print_function, absolute_import, division

from KratosMultiphysics.SolidMechanicsApplication import *
from KratosMultiphysics.MeshingApplication import *
from KratosMultiphysics.PfemFluidDynamicsApplication import *
from KratosMultiphysics.DEMApplication import *

# Application dependent names and paths
from KratosMultiphysics import _ImportApplication
from KratosFemToDemApplication import *
application = KratosFemToDemApplication()
application_name = "KratosFemToDemApplication"

_ImportApplication(application, application_name)