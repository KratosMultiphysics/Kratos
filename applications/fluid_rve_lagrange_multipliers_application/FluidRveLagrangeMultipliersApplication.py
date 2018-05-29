from KratosFluidRveLagrangeMultipliersApplication import *
application = KratosFluidRveLagrangeMultipliersApplication()
application_name = "KratosFluidRveLagrangeMultipliersApplication"
application_folder = "fluid_rve_lagrange_multipliers_application"

# The following lines are common for all applications
from . import application_importer
import inspect
caller = inspect.stack()[1] # Information about the file that imported this, to check for unexpected imports
application_importer.ImportApplication(application,application_name,application_folder,caller)
