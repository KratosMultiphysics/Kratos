# Application dependent names and paths
from KratosMultiphysics import _ImportApplication, python_registry_utilities
from KratosSPHApplication import *

application = KratosSPHApplication()
application_name = "KratosSPHApplication"

_ImportApplication(application, application_name)

from . import python_registry_lists
python_registry_utilities.RegisterModelersList("KratosMultiphysics.SPHApplication", python_registry_lists)
python_registry_utilities.RegisterOperationsList("KratosMultiphysics.SPHApplication", python_registry_lists)
python_registry_utilities.RegisterProcessesList("KratosMultiphysics.SPHApplication", python_registry_lists)