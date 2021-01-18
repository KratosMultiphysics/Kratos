# Application dependent names and paths
from KratosMultiphysics import _ImportApplication
from KratosStructuralMechanicsApplication import *
application = KratosStructuralMechanicsApplication()
application_name = "KratosStructuralMechanicsApplication"

_ImportApplication(application, application_name)
