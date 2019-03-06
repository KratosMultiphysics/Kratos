# ==============================================================================
#  KratosShapeOptimizationApplication
#
#  License:         BSD License
#                   license: ShapeOptimizationApplication/license.txt
#
#  Main authors:    Baumgaertner Daniel, https://github.com/dbaumgaertner
#
# ==============================================================================

# ------------------------------------------------------------------------------
# Imports
# ------------------------------------------------------------------------------
from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
from KratosShapeOptimizationApplication import *
application = KratosShapeOptimizationApplication()
application_name = "KratosShapeOptimizationApplication"
application_folder = "ShapeOptimizationApplication"

# The following lines are common for all applications
from .. import application_importer
import inspect
caller = inspect.stack()[1] # Information about the file that imported this, to check for unexpected imports
application_importer.ImportApplication(application,application_name,application_folder,caller,__path__)
