from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

import KratosMultiphysics
from KratosConvectionDiffusionApplication import *
application = KratosConvectionDiffusionApplication()
application_name = "KratosConvectionDiffusionApplication"
application_folder = "ConvectionDiffusionApplication"

KratosMultiphysics._ImportApplicationAsModule(application, application_name, application_folder, __path__)