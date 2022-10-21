# Application dependent names and paths
from KratosMultiphysics import _ImportApplication
from KratosPFEM2Application import *
application = KratosPFEM2Application()
application_name = "KratosPFEM2Application"

_ImportApplication(application,application_name)
