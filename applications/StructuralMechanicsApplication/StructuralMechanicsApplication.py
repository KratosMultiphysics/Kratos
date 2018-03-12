from KratosStructuralMechanicsApplication import *
application = KratosStructuralMechanicsApplication()
application_name = "KratosStructuralMechanicsApplication"
application_folder = "StructuralMechanicsApplication"

# The following lines are common for all applications
from . import application_importer
import inspect
caller = inspect.stack()[1] # Information about the file that imported this, to check for unexpected imports
application_importer.ImportApplication(application,application_name,application_folder,caller)

# Add custom_response path to system path, so that the python scripts can be found
import os
import sys
import KratosMultiphysics
Globals = KratosMultiphysics.KratosGlobals
applications_root = Globals.ApplicationsRoot
application_path = os.path.join(applications_root, application_folder)
python_path = os.path.join(application_path, 'custom_response_functions')
sys.path.append(python_path)