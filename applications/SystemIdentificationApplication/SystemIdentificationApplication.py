# Application dependent names and paths
from KratosMultiphysics import _ImportApplication, python_registry_utilities
from KratosSystemIdentificationApplication import *

application = KratosSystemIdentificationApplication()
application_name = "KratosSystemIdentificationApplication"

_ImportApplication(application, application_name)

from . import python_registry_lists
python_registry_utilities.RegisterModelersList("KratosMultiphysics.SystemIdentificationApplication", python_registry_lists)
python_registry_utilities.RegisterOperationsList("KratosMultiphysics.SystemIdentificationApplication", python_registry_lists)
python_registry_utilities.RegisterProcessesList("KratosMultiphysics.SystemIdentificationApplication", python_registry_lists)