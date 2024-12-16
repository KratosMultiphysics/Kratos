
import KratosMultiphysics.PfemFluidDynamicsApplication
import KratosMultiphysics.DEMApplication
import KratosMultiphysics.StructuralMechanicsApplication
import KratosMultiphysics.ConstitutiveLawsApplication

# Application dependent names and paths
from KratosMultiphysics import _ImportApplication
from KratosFemToDemApplication import *
application = KratosFemToDemApplication()
application_name = "KratosFemToDemApplication"

_ImportApplication(application, application_name)