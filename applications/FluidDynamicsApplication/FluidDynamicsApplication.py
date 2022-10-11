from KratosFluidDynamicsApplication import *

from KratosMultiphysics import _ImportApplication
application = KratosFluidDynamicsApplication()
application_name = "KratosFluidDynamicsApplication"

_ImportApplication(application, application_name)
