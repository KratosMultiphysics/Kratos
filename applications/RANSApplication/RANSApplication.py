# makes KratosMultiphysics backward compatible with python 2.6 and 2.7
from __future__ import print_function, absolute_import, division

# Application dependent names and paths
from KratosMultiphysics.kratos_utilities import CheckIfApplicationsAvailable
from KratosMultiphysics import IsDistributedRun
if (IsDistributedRun() and CheckIfApplicationsAvailable("TrilinosApplication")):
    from KratosMultiphysics.TrilinosApplication import *
from KratosRANSApplication import *

from KratosMultiphysics import _ImportApplicationAsModule
application = KratosRANSApplication()
application_name = "KratosRANSApplication"
application_folder = "RANSApplication"

_ImportApplicationAsModule(application, application_name, application_folder, __path__)
