# makes KratosMultiphysics backward compatible with python 2.6 and 2.7
from __future__ import print_function, absolute_import, division

# Application dependent names and paths
from KratosIGAStructuralMechanicsApplication import *
application = KratosIGAStructuralMechanicsApplication()
application_name = "KratosIGAStructuralMechanicsApplication"
application_folder = "IGAStructuralMechanicsApplication"

# The following lines are common for all applications
from . import application_importer
import inspect
caller = inspect.stack()[1] # Information about the file that imported this, to check for unexpected imports
application_importer.ImportApplication(application,application_name,application_folder,caller)
