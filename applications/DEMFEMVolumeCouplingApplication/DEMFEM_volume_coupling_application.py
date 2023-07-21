# Application dependent names and paths
from KratosMultiphysics import _ImportApplication
from KratosDEMFEMVolumeCouplingApplication import *
application = KratosDEMFEMVolumeCouplingApplication()
application_name = "KratosDEMFEMVolumeCouplingApplication"

_ImportApplication(application, application_name)