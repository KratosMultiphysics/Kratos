
# Application dependent names and paths
from KratosMultiphysics import _ImportApplication
from KratosConstitutiveModelsApplication import *
application = KratosConstitutiveModelsApplication()
application_name = "KratosConstitutiveModelsApplication"

_ImportApplication(application, application_name)