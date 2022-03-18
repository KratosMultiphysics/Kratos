# Import the applications the ThermalDEMApplication depends on
import KratosMultiphysics.DEMApplication

# Application dependent names and paths
from KratosMultiphysics import _ImportApplication
from KratosThermalDEMApplication import *
application = KratosThermalDEMApplication()
application_name = "KratosThermalDEMApplication"

_ImportApplication(application, application_name)
