# makes KratosMultiphysics backward compatible with python 2.6 and 2.7
from __future__ import print_function, absolute_import, division

# Application dependent names and paths
from KratosMultiphysics import _ImportApplication

# why is this import necessary here ?
import KratosMultiphysics.FluidDynamicsApplication
# import KratosMultiphysics.RANSApplication

from KratosChimeraApplication import *

application = KratosChimeraApplication()
application_name = "KratosChimeraApplication"

_ImportApplication(application, application_name)
