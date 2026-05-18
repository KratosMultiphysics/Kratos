# Application dependent names and paths
import KratosMultiphysics
from KratosMultiphysics import _ImportApplication, python_registry_utilities
import KratosMultiphysics.StructuralMechanicsApplication
import KratosGeoMechanicsApplication as _KratosGeoMechanicsApplication

# Explicitly export required bindings by name.
_required_bindings = (
    "KratosGeoMechanicsApplication",
    "ProcessUtilities",
    "CustomWorkflowFactory",
    "UMAT_PARAMETERS",
)
for _name in _required_bindings:
    globals()[_name] = getattr(_KratosGeoMechanicsApplication, _name)

KratosGeoMechanicsApplication = _KratosGeoMechanicsApplication.KratosGeoMechanicsApplication

# Re-export all public bindings from the extension module without wildcard imports.
for _name in dir(_KratosGeoMechanicsApplication):
    if not _name.startswith("_"):
        globals()[_name] = getattr(_KratosGeoMechanicsApplication, _name)

application = KratosGeoMechanicsApplication()
application_name = "KratosGeoMechanicsApplication"

_ImportApplication(application, application_name)

if not KratosMultiphysics.Registry.HasItem("Stages.KratosMultiphysics.GeoMechanicsApplication.GeoMechanicsAnalysis"):
    from . import python_registry_lists
    python_registry_utilities.RegisterAll("KratosMultiphysics.GeoMechanicsApplication", python_registry_lists)