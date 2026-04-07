
# Application dependent names and paths
import KratosMultiphysics
from KratosMultiphysics import _ImportApplication
from KratosViscosityModulatorApplication import *
application = KratosViscosityModulatorApplication()
application_name = "KratosViscosityModulatorApplication"

_ImportApplication(application, application_name)
