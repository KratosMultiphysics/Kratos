#
#   KRATOS  _____________
#          /  _/ ____/   |
#          / // / __/ /| |
#        _/ // /_/ / ___ |
#       /___/\____/_/  |_| Application
#
#   Main authors:   Thomas Oberbichler
#

# makes KratosMultiphysics backward compatible with python 2.6 and 2.7
from __future__ import print_function, absolute_import, division

# Application dependent names and paths
from KratosMultiphysics import _ImportApplication
from KratosIgaApplication import *
application = KratosIgaApplication()
application_name = "KratosIgaApplication"

_ImportApplication(application, application_name)
