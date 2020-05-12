# makes KratosMultiphysics backward compatible with python 2.6 and 2.7
from __future__ import print_function, absolute_import, division

import sys, os

# Application dependent names and paths
from KratosMultiphysics import _ImportApplication
from KratosMultilevelMonteCarloApplication import *
application = KratosMultilevelMonteCarloApplication()
application_name = "KratosMultilevelMonteCarloApplication"

_ImportApplication(application, application_name)

sys.path.append(os.path.join(os.path.dirname(sys.path[0]),'XMC'))
