
# Application dependent names and paths
from KratosMultiphysics import _ImportApplication
from KratosRomApplication import *
application = KratosRomApplication()
application_name = "KratosRomApplication"

_ImportApplication(application, application_name)
