# Application dependent names and paths
from KratosMultiphysics import _ImportApplication
from KratosConstitutiveLawsApplication import *
application = KratosConstitutiveLawsApplication()
application_name = "KratosConstitutiveLawsApplication"

_ImportApplication(application, application_name)
