# makes KratosMultiphysics backward compatible with python 2.6 and 2.7
from __future__ import print_function, absolute_import, division

# Application dependent names and paths
import KratosMultiphysics as KM
from KratosCoSimulationApplication import *
application = KratosCoSimulationApplication()
application_name = "KratosCoSimulationApplication"
application_folder = "CoSimulationApplication"

KM._ImportApplicationAsModule(application, application_name, application_folder, __path__)
