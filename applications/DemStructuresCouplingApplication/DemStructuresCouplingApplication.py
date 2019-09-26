# makes KratosMultiphysics backward compatible with python 2.6 and 2.7
from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Application dependent names and paths
import KratosMultiphysics as KM
from KratosDemStructuresCouplingApplication import *
application = KratosDemStructuresCouplingApplication()
application_name = "KratosDemStructuresCouplingApplication"
application_folder = "DemStructuresCouplingApplication"

KM._ImportApplicationAsModule(application, application_name, application_folder, __path__)