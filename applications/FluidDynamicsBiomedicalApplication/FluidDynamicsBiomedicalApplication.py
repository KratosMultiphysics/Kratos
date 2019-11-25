from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
from KratosFluidDynamicsBiomedicalApplication import *

from KratosMultiphysics import _ImportApplicationAsModule
application = KratosFluidDynamicsBiomedicalApplication()
application_name = "KratosFluidDynamicsBiomedicalApplication"
application_folder = "FluidDynamicsBiomedicalApplication"

_ImportApplicationAsModule(application, application_name, application_folder, __path__)
