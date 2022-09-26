
# Application dependent names and paths
from KratosMultiphysics import _ImportApplication
from KratosPoromechanicsApplication import *
application = KratosPoromechanicsApplication()
application_name = "KratosPoromechanicsApplication"

_ImportApplication(application, application_name)