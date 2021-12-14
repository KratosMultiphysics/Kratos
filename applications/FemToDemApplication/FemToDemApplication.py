
# Application dependent names and paths
from KratosMultiphysics import _ImportApplication
from KratosFemToDemApplication import *
application = KratosFemToDemApplication()
application_name = "KratosFemToDemApplication"

_ImportApplication(application, application_name)