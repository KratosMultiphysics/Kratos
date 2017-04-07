# ==============================================================================
#  TopologyOptimizationApplication
#
#  License:         BSD License
#                   license: TopologyOptimizationApplication/license.txt
#
#  Main authors:    Baumgärtner Daniel, https://github.com/dbaumgaertner
#
# ==============================================================================

# ------------------------------------------------------------------------------
# Imports
# ------------------------------------------------------------------------------

from KratosTopologyOptimizationApplication import *
application = KratosTopologyOptimizationApplication()
application_name = "KratosTopologyOptimizationApplication"
application_folder = "TopologyOptimizationApplication"

# The following lines are common for all applications
from . import application_importer
import inspect
caller = inspect.stack()[1] # Information about the file that imported this, to check for unexpected imports
application_importer.ImportApplication(application,application_name,application_folder,caller)
