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
from KratosMultiphysics import _ImportApplication
from KratosShapeOptimizationApplication import *
application = KratosShapeOptimizationApplication()
application_name = "KratosShapeOptimizationApplication"

_ImportApplication(application, application_name)