# Application dependent names and paths
from KratosMultiphysics import _ImportApplication
from KratosStatisticsApplication import *
application = KratosStatisticsApplication()
application_name = "KratosStatisticsApplication"

_ImportApplication(application, application_name)
