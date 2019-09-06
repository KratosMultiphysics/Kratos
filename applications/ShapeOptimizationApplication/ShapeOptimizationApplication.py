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
import KratosMultiphysics as KM
from KratosShapeOptimizationApplication import *
application = KratosShapeOptimizationApplication()
application_name = "KratosShapeOptimizationApplication"
application_folder = "ShapeOptimizationApplication"

KM._ImportApplicationAsModule(application, application_name, application_folder, __path__)