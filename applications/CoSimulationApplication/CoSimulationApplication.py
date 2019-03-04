# makes KratosMultiphysics backward compatible with python 2.6 and 2.7
from __future__ import print_function, absolute_import, division
import os
import sys

# Application dependent names and paths
import KratosMultiphysics as KM
from KratosCoSimulationApplication import *
application = KratosCoSimulationApplication()
application_name = "KratosCoSimulationApplication"
application_folder = "CoSimulationApplication"
folder_path  = __path__

from .. import application_importer
import inspect
caller = inspect.stack()[1] # Information about the file that imported this, to check for unexpected imports
application_importer.ImportApplication(application,application_name,application_folder,caller, __path__)

## adding pyKratos to path
python_path = os.path.join(folder_path[0], 'custom_data_structure')
sys.path.append(python_path)
