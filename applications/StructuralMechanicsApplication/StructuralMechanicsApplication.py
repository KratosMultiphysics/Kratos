# Application dependent names and paths
from KratosMultiphysics import _ImportApplication
from KratosMultiphysics.kratos_utilities import CheckIfApplicationsAvailable
from KratosStructuralMechanicsApplication import *
application = KratosStructuralMechanicsApplication()
application_name = "KratosStructuralMechanicsApplication"

_ImportApplication(application, application_name)

if CheckIfApplicationsAvailable("ConstitutiveLawsApplication"):
    # if available import the advanced constitutive laws
    import KratosMultiphysics.ConstitutiveLawsApplication
