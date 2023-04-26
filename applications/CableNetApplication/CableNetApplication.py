
# Application dependent names and paths
from KratosMultiphysics import _ImportApplication
import KratosMultiphysics.StructuralMechanicsApplication
from KratosCableNetApplication import *
application = KratosCableNetApplication()
application_name = "KratosCableNetApplication"

_ImportApplication(application, application_name)
