
# Application dependent names and paths
from KratosMultiphysics import _ImportApplication
import KratosMultiphysics.ConstitutiveModelsApplication
from KratosUmatApplication import *
application = KratosUmatApplication()
application_name = "KratosUmatApplication"

_ImportApplication(application, application_name)
