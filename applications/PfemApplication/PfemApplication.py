# Applications requiered
from KratosMultiphysics.DelaunayMeshingApplication import *
from KratosMultiphysics.SolidMechanicsApplication import *

# Application dependent names and paths
from KratosMultiphysics import _ImportApplication
from KratosPfemApplication import *
application = KratosPfemApplication()
application_name = "KratosPfemApplication"

_ImportApplication(application, application_name)
