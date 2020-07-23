from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
from KratosSwimmingDEMApplication import *

from KratosMultiphysics import _ImportApplication
application = KratosSwimmingDEMApplication()
application_name = "KratosSwimmingDEMApplication"

_ImportApplication(application, application_name)
