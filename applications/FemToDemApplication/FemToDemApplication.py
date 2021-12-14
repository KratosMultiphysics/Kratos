
import KratosMultiphysics.MeshingApplication
import KratosMultiphysics.PfemFluidDynamicsApplication
import KratosMultiphysics.DEMApplication
import KratosMultiphysics.StructuralMechanicsApplication

# Application dependent names and paths
from KratosMultiphysics import _ImportApplication
from KratosFemToDemApplication import *
application = KratosFemToDemApplication()
application_name = "KratosFemToDemApplication"

_ImportApplication(application, application_name)