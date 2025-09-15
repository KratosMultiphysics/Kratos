
# Application dependent names and paths
from KratosMultiphysics import _ImportApplication
from KratosConvectionDiffusionApplication import *
application = KratosConvectionDiffusionApplication()
application_name = "KratosConvectionDiffusionApplication"

_ImportApplication(application, application_name)
