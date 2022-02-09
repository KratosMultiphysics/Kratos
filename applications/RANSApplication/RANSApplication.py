# Application dependent names and paths
from KratosMultiphysics.kratos_utilities import CheckIfApplicationsAvailable
from KratosMultiphysics import IsDistributedRun
if (IsDistributedRun() and CheckIfApplicationsAvailable("TrilinosApplication")):
    import KratosMultiphysics.TrilinosApplication
elif (IsDistributedRun()):
    raise Exception("Distributed run requires TrilinosApplication")

import KratosMultiphysics.FluidDynamicsApplication
from KratosRANSApplication import *

from KratosMultiphysics import _ImportApplication
application = KratosRANSApplication()
application_name = "KratosRANSApplication"

_ImportApplication(application, application_name)
