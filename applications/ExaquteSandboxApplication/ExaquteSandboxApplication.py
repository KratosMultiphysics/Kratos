# makes KratosMultiphysics backward compatible with python 2.6 and 2.7
from __future__ import print_function, absolute_import, division

from KratosMultiphysics import _ImportApplication
from KratosExaquteSandboxApplication import *
application = KratosExaquteSandboxApplication()
application_name = "KratosExaquteSandboxApplication"

_ImportApplication(application, application_name)