from KratosExternalSolversApplication import *
application = KratosExternalSolversApplication()
application_name = "KratosExternalSolversApplication"
application_folder = "ExternalSolversApplication"

# The following lines are common for all applications
import application_importer
import inspect
caller = inspect.stack()[1] # Information about the file that imported this, to check for unexpected imports
application_importer.ImportApplication(application,application_name,application_folder,caller)
