# Application dependent names and paths
from KratosMultiphysics import _ImportApplication, python_registry_utilities
from KratosMedApplication import *

application = KratosMedApplication()
application_name = "KratosMedApplication"

_ImportApplication(application, application_name)

from . import python_registry_lists

python_registry_utilities.RegisterAll("KratosMultiphysics.MedApplication", python_registry_lists)
