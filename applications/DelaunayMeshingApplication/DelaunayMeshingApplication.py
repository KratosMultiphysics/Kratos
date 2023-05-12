
from KratosMultiphysics import _ImportApplication

from KratosDelaunayMeshingApplication import *
application = KratosDelaunayMeshingApplication()
application_name = "KratosDelaunayMeshingApplication"

_ImportApplication(application, application_name)
