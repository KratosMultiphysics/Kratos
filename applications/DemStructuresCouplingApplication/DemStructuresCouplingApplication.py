# makes KratosMultiphysics backward compatible with python 2.6 and 2.7
from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Application dependent names and paths
import KratosMultiphysics as KM
from KratosDemStructuresCouplingApplication import *
application = KratosDemStructuresCouplingApplication()
application_name = "KratosDemStructuresCouplingApplication"

KM._ImportApplication(application, application_name)