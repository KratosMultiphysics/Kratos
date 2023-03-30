
# Application dependent names and paths
from KratosMultiphysics import _ImportApplication
import KratosMultiphysics.StructuralMechanicsApplication
import KratosMultiphysics.PoromechanicsApplication
from KratosDamApplication import *
application = KratosDamApplication()
application_name = "KratosDamApplication"

_ImportApplication(application,application_name)
