# makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Application dependent names and paths
from KratosMultiphysics import _ImportApplication
from KratosConstitutiveModelsApplication import *
application = KratosConstitutiveModelsApplication()
application_name = "KratosConstitutiveModelsApplication"

_ImportApplication(application, application_name)