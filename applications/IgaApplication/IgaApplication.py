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
from KratosIgaApplication import *
application = KratosIgaApplication()
application_name = "KratosIgaApplication"
application_folder = "IgaApplication"

# The following lines are common for all applications
from . import application_importer
import inspect
# Information about the file that imported this, to check for unexpected imports
caller = inspect.stack()[1]
application_importer.ImportApplication(application, application_name,
    application_folder, caller)
