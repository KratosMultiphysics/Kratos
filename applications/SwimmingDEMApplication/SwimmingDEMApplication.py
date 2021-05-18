from KratosMultiphysics import _ImportApplication

# import the applications the SwimmingDEMApplication depends on
import KratosMultiphysics.FluidDynamicsApplication
import KratosMultiphysics.DEMApplication

from KratosSwimmingDEMApplication import *
application = KratosSwimmingDEMApplication()
application_name = "KratosSwimmingDEMApplication"

_ImportApplication(application, application_name)
