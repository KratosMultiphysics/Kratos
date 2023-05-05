import os

# Application dependent names and paths
from KratosMultiphysics import _ImportApplication
import KratosMultiphysics.StructuralMechanicsApplication
from KratosGeoMechanicsApplication import *
application = KratosGeoMechanicsApplication()
application_name = "KratosGeoMechanicsApplication"

geo_experimental = os.getenv('GEO_EXPERIMENTAL')
if geo_experimental:
    application.SetExperimental()

_ImportApplication(application, application_name)
