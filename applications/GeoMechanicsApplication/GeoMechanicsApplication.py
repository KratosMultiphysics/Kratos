# Application dependent names and paths
from KratosMultiphysics import _ImportApplication, python_registry_utilities
import KratosMultiphysics.StructuralMechanicsApplication
import KratosGeoMechanicsApplication as _KratosGeoMechanicsApplication

# Keep module-level exports equivalent to the previous wildcard import.
for _name in dir(_KratosGeoMechanicsApplication):
    if not _name.startswith("_"):
        globals()[_name] = getattr(_KratosGeoMechanicsApplication, _name)

application = KratosGeoMechanicsApplication()
application_name = "KratosGeoMechanicsApplication"

_ImportApplication(application, application_name)

if not KratosMultiphysics.Registry.HasItem("Stages.KratosMultiphysics.GeoMechanicsApplication.GeoMechanicsAnalysis"):
    from . import python_registry_lists
    python_registry_utilities.RegisterAll("KratosMultiphysics.GeoMechanicsApplication", python_registry_lists)