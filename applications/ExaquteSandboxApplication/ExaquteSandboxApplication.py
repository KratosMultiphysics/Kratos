# makes KratosMultiphysics backward compatible with python 2.6 and 2.7
from __future__ import print_function, absolute_import, division

import KratosMultiphysics as KM
from KratosExaquteSandboxApplication import *
application = KratosExaquteSandboxApplication()
application_name = "KratosExaquteSandboxApplication"
application_folder = "ExaquteSandboxApplication"

KM._ImportApplicationAsModule(application, application_name, application_folder, __path__)