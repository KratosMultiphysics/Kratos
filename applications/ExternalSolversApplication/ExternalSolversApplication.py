# makes KratosMultiphysics backward compatible with python 2.6 and 2.7
from __future__ import print_function, absolute_import, division

# Application dependent names and paths
import KratosMultiphysics as KM
from KratosExternalSolversApplication import *
application = KratosExternalSolversApplication()
application_name = "ExternalSolversApplication"

KM._ImportApplication(application, application_name)
