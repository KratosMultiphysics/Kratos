from KratosMultiphysics import _ImportApplication
from KratosWindEngineeringApplication import *

application = KratosWindEngineeringApplication()
application_name = "KratosWindEngineeringApplication"

_ImportApplication(application, application_name)