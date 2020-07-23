# makes KratosMultiphysics backward compatible with python 2.6 and 2.7
from __future__ import print_function, absolute_import, division

# Applications requiered
from KratosMultiphysics.ConstitutiveModelsApplication import *

# Application dependent names and paths
from KratosMultiphysics import _ImportApplication
from KratosSolidMechanicsApplication import *
application = KratosSolidMechanicsApplication()
application_name = "KratosSolidMechanicsApplication"

_ImportApplication(application, application_name)