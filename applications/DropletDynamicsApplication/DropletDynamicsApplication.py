# Application dependent names and paths
from KratosMultiphysics import _ImportApplication
from KratosDropletDynamicsApplication import *
application = KratosDropletDynamicsApplication()
application_name = "KratosDropletDynamicsApplication"

_ImportApplication(application, application_name)
