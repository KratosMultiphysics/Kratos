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
import KratosMultiphysics as KM
from KratosIgaApplication import *
application = KratosIgaApplication()
application_name = "KratosIgaApplication"
application_folder = "IgaApplication"

KM._ImportApplicationAsModule(application, application_name, application_folder, __path__)
