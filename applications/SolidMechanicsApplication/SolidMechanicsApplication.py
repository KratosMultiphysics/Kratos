# Applications required
from KratosMultiphysics.ConstitutiveModelsApplication import *

# Application dependent names and paths
from KratosMultiphysics import _ImportApplication
from KratosSolidMechanicsApplication import *
application = KratosSolidMechanicsApplication()
application_name = "KratosSolidMechanicsApplication"

_ImportApplication(application, application_name)
