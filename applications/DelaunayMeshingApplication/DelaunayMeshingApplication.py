from __future__ import print_function, absolute_import, division

import KratosMultiphysics as KM

from KratosDelaunayMeshingApplication import *
application = KratosDelaunayMeshingApplication()
application_name = "KratosDelaunayMeshingApplication"
application_folder = "DelaunayMeshingApplication"

KM._ImportApplicationAsModule(application, application_name, application_folder, __path__)
