from KratosMultiphysics import _ImportApplication
from KratosKaHIPApplication import *
application = KratosKaHIPApplication()
application_name = "KratosKaHIPApplication"

_ImportApplication(application, application_name)
