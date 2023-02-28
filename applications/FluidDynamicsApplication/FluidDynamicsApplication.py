from KratosFluidDynamicsApplication import *

from KratosMultiphysics import _ImportApplication, python_registry_utilities

application = KratosFluidDynamicsApplication()
application_name = "KratosFluidDynamicsApplication"

_ImportApplication(application, application_name)

from . import python_registry_lists
python_registry_utilities.RegisterAll("KratosMultiphysics.FluidDynamicsApplication", python_registry_lists)
