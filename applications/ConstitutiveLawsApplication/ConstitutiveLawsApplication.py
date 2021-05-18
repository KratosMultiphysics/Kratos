# Application dependent names and paths
from KratosMultiphysics import _ImportApplication
import KratosMultiphysics.StructuralMechanicsApplication # temp until dependencies are moved to Core
from KratosConstitutiveLawsApplication import *
application = KratosConstitutiveLawsApplication()
application_name = "KratosConstitutiveLawsApplication"

_ImportApplication(application, application_name)
