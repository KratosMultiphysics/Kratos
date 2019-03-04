# Applications requiered
from KratosMultiphysics.ContactMechanicsApplication import *
from KratosMultiphysics.SolidMechanicsApplication import *

from KratosPfemSolidMechanicsApplication import *
application = KratosPfemSolidMechanicsApplication()
application_name = "KratosPfemSolidMechanicsApplication"
application_folder = "PfemSolidMechanicsApplication"

# The following lines are common for all applications
from .. import application_importer
import inspect
caller = inspect.stack()[1] # Information about the file that imported this, to check for unexpected imports
application_importer.ImportApplication(application,application_name,application_folder,caller, __path__)
