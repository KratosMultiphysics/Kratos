# Application dependent names and paths
from KratosMultiphysics import _ImportApplication, python_registry_utilities
import KratosMultiphysics.StructuralMechanicsApplication
from KratosGeoMechanicsApplication import *
application = KratosGeoMechanicsApplication()
application_name = "KratosGeoMechanicsApplication"

_ImportApplication(application, application_name)

if not KratosMultiphysics.Registry.HasItem("Stages.KratosMultiphysics.GeoMechanicsApplication.GeoMechanicsAnalysis"):
    from . import python_registry_lists
    python_registry_utilities.RegisterAll("KratosMultiphysics.GeoMechanicsApplication", python_registry_lists)