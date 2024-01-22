from KratosMultiphysics import _ImportApplication
from KratosParticleMechanicsApplication import *
application = KratosParticleMechanicsApplication()
application_name = "KratosParticleMechanicsApplication"

_ImportApplication(application, application_name)
