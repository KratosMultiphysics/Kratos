from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
from KratosFluidDynamicsApplication import *

from KratosMultiphysics import _ImportApplication
application = KratosFluidDynamicsApplication()
application_name = "KratosFluidDynamicsApplication"

_ImportApplication(application, application_name)
