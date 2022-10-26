# Application dependent names and paths
from KratosMultiphysics import _ImportApplication
from KratosMedApplication import *
application = KratosMedApplication()
application_name = "KratosMedApplication"

_ImportApplication(application, application_name)
