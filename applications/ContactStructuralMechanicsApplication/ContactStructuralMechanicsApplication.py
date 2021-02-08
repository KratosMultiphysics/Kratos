# Application dependent names and paths
from KratosMultiphysics import _ImportApplication
from KratosContactStructuralMechanicsApplication import *
application = KratosContactStructuralMechanicsApplication()
application_name = "KratosContactStructuralMechanicsApplication"

_ImportApplication(application, application_name)
