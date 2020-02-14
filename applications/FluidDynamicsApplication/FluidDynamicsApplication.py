from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
from KratosFluidDynamicsApplication import *

from KratosMultiphysics import _ImportApplicationAsModule
application = KratosFluidDynamicsApplication()
application_name = "KratosFluidDynamicsApplication"
application_folder = "FluidDynamicsApplication"

_ImportApplicationAsModule(application, application_name, application_folder, __path__)
