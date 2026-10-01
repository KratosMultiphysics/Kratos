# Application dependent names and paths
from KratosMultiphysics import _ImportApplication, python_registry_utilities
from KratosMeshioPlusPlusApplication import *

application = KratosMeshioPlusPlusApplication()
application_name = "KratosMeshioPlusPlusApplication"

_ImportApplication(application, application_name)

from . import python_registry_lists

python_registry_utilities.RegisterAll("KratosMultiphysics.MeshioPlusPlusApplication", python_registry_lists)
