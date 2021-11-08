# Applications requiered
from KratosMultiphysics.ContactMechanicsApplication import *
from KratosMultiphysics.PfemApplication import *

# Application dependent names and paths
from KratosMultiphysics import _ImportApplication
from KratosPfemSolidMechanicsApplication import *
application = KratosPfemSolidMechanicsApplication()
application_name = "KratosPfemSolidMechanicsApplication"

_ImportApplication(application, application_name)
