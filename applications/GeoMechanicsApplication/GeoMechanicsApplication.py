# Application dependent names and paths
from KratosMultiphysics import _ImportApplication
import KratosMultiphysics.StructuralMechanicsApplication
from KratosGeoMechanicsApplication import *
application = KratosGeoMechanicsApplication()
application_name = "KratosGeoMechanicsApplication"

_ImportApplication(application, application_name)
