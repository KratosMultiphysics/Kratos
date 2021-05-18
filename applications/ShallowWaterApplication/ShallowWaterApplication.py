# Application dependent names and paths
from KratosMultiphysics import _ImportApplication
from KratosShallowWaterApplication import *
application = KratosShallowWaterApplication()
application_name = "KratosShallowWaterApplication"

_ImportApplication(application, application_name)
