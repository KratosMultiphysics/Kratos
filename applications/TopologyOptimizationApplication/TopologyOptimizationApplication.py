# Application dependent names and paths
from KratosMultiphysics import _ImportApplication
from KratosTopologyOptimizationApplication import *
application = KratosTopologyOptimizationApplication()
application_name = "KratosTopologyOptimizationApplication"

_ImportApplication(application, application_name)
