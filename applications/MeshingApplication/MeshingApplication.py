# Application dependent names and paths
from KratosMultiphysics import _ImportApplication
from KratosMeshingApplication import *
application = KratosMeshingApplication()
application_name = "KratosMeshingApplication"

_ImportApplication(application, application_name)
