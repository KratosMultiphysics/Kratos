# makes KratosMultiphysics backward compatible with python 2.6 and 2.7
from __future__ import print_function, absolute_import, division

# Application dependent names and paths
from KratosMultiphysics import _ImportApplication
import KratosMultiphysics.FluidDynamicsApplication
from KratosChimeraApplication import *
application = KratosChimeraApplication()
application_name = "KratosChimeraApplication"

_ImportApplication(application, application_name)
