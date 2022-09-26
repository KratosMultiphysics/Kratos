
# Application dependent names and paths
from KratosMultiphysics import _ImportApplication
from KratosFreeSurfaceApplication import *
application = KratosFreeSurfaceApplication()
application_name = "KratosFreeSurfaceApplication"

_ImportApplication(application, application_name)
