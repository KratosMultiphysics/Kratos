from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
from KratosSwimmingDEMApplication import *

from KratosMultiphysics import _ImportApplicationAsModule
application = KratosSwimmingDEMApplication()
application_name = "KratosSwimmingDEMApplication"
application_folder = "SwimmingDEMApplication"

_ImportApplicationAsModule(application, application_name, application_folder, __path__)
