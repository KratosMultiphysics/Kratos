# makes KratosMultiphysics backward compatible with python 2.6 and 2.7
from __future__ import print_function, absolute_import, division

# Application dependent names and paths
from KratosRANSModellingApplication import *

from KratosMultiphysics import _ImportApplicationAsModule
application = KratosRANSModellingApplication()
application_name = "KratosRANSModellingApplication"
application_folder = "RANSModellingApplication"

_ImportApplicationAsModule(application, application_name, application_folder, __path__)
