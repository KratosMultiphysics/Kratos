# Applications requiered
from KratosMultiphysics.SolversApplication import *
from KratosMultiphysics.ConstitutiveModelsApplication import *

from KratosSolidMechanicsApplication import *
application = KratosSolidMechanicsApplication()
application_name = "KratosSolidMechanicsApplication"
application_folder = "SolidMechanicsApplication"

# The following lines are common for all applications
from .. import application_importer
import inspect
caller = inspect.stack()[1] # Information about the file that imported this, to check for unexpected imports
application_importer.ImportApplication(application,application_name,application_folder,caller, __path__)
