#
#   KRATOS  _____________
#          /  _/ ____/   |
#          / // / __/ /| |
#        _/ // /_/ / ___ |
#       /___/\____/_/  |_| Application
#
#   Main authors:   Thomas Oberbichler
#

# Application dependent names and paths
from KratosMultiphysics import _ImportApplication
from KratosIgaApplication import *
application = KratosIgaApplication()
application_name = "KratosIgaApplication"

_ImportApplication(application, application_name)
