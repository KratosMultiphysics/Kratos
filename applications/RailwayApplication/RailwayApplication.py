# Application dependent names and paths
from KratosMultiphysics import _ImportApplication, python_registry_utilities
import KratosMultiphysics.StructuralMechanicsApplication
import KratosMultiphysics.GeomechanicsApplication
from KratosRailwayApplication import *
application = KratosRailwayApplication()
application_name = "KratosRailwayApplication"

_ImportApplication(application, application_name)
