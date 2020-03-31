from __future__ import print_function, absolute_import, division

import KratosMultiphysics as KM

from KratosDelaunayMeshingApplication import *
application = KratosDelaunayMeshingApplication()
application_name = "KratosDelaunayMeshingApplication"

KM._ImportApplication(application, application_name)
