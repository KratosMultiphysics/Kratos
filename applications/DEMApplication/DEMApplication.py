from KratosDEMApplication import *

from KratosMultiphysics import _ImportApplication
application = KratosDEMApplication()
application_name = "KratosDEMApplication"

_ImportApplication(application, application_name)
