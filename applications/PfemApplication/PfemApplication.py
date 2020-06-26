# Applications requiered
from KratosMultiphysics.DelaunayMeshingApplication import *
from KratosMultiphysics.SolidMechanicsApplication import *

from KratosPfemApplication import *
application = KratosPfemApplication()
application_name = "KratosPfemApplication"
application_folder = "PfemApplication"

# The following lines are common for all applications
from .. import application_importer
import inspect
caller = inspect.stack()[1] # Information about the file that imported this, to check for unexpected imports
application_importer.ImportApplication(application,application_name,application_folder,caller, __path__)
