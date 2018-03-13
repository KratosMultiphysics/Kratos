# 
#   KRATOS _______
#         / ____(_)___ ____  ____
#        / __/ / / __ `/ _ \/ __ \
#       / /___/ / /_/ /  __/ / / /
#      /_____/_/\__, /\___/_/ /_/ SolversApplication
#              /____/
# 
#   Author: Thomas Oberbichler
#

from KratosEigenSolversApplication import *
application = KratosEigenSolversApplication()
application_name = "KratosEigenSolversApplication"
application_folder = "eigen_solvers_application"

# The following lines are common for all applications
from . import application_importer
import inspect
caller = inspect.stack()[1] # Information about the file that imported this, to check for unexpected imports
application_importer.ImportApplication(application,application_name,application_folder,caller)
