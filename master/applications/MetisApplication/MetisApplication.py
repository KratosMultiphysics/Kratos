from KratosMultiphysics import _ImportApplication
from KratosMetisApplication import *
application = KratosMetisApplication()
application_name = "KratosMetisApplication"

_ImportApplication(application, application_name)
