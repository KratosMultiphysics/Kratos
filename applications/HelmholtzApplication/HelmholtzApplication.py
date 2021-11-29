# Application dependent names and paths
from KratosMultiphysics import _ImportApplication
from KratosHelmholtzApplication import *
application = KratosHelmholtzApplication()
application_name = "KratosHelmholtzApplication"

_ImportApplication(application, application_name)
