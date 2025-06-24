# Applications required
from KratosMultiphysics.DelaunayMeshingApplication import *

# Application dependent names and paths
from KratosMultiphysics import _ImportApplication
import KratosMultiphysics.DelaunayMeshingApplication
from KratosContactMechanicsApplication import *
application = KratosContactMechanicsApplication()
application_name = "KratosContactMechanicsApplication"

_ImportApplication(application, application_name)
