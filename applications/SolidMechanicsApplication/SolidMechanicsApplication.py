# makes KratosMultiphysics backward compatible with python 2.6 and 2.7
from __future__ import print_function, absolute_import, division

# Applications requiered
from KratosMultiphysics.ConstitutiveModelsApplication import *

# Application dependent names and paths
import KratosMultiphysics as KM
from KratosSolidMechanicsApplication import *
application = KratosSolidMechanicsApplication()
application_name = "KratosSolidMechanicsApplication"
application_folder = "SolidMechanicsApplication"

KM._ImportApplicationAsModule(application, application_name, application_folder, __path__)