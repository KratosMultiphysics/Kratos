# makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Application dependent names and paths
from KratosMultiphysics import _ImportApplication
from KratosRomApplication import *
application = KratosRomApplication()
application_name = "KratosRomApplication"

_ImportApplication(application, application_name)
