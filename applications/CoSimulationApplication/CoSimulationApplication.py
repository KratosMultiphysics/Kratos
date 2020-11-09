# Application dependent names and paths
from KratosMultiphysics import _ImportApplication
from KratosCoSimulationApplication import *
application = KratosCoSimulationApplication()
application_name = "KratosCoSimulationApplication"

_ImportApplication(application, application_name)
