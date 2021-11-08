from __future__ import print_function, absolute_import, division

from KratosMultiphysics import _ImportApplication

from KratosDelaunayMeshingApplication import *
application = KratosDelaunayMeshingApplication()
application_name = "KratosDelaunayMeshingApplication"

_ImportApplication(application, application_name)
