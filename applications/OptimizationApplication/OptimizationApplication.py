# ==============================================================================
#  KratosOptimizationApplication
#
#  License:         BSD License
#                   license: OptimizationApplication/license.txt
#
#  Main authors:    Reza Najian Asl, https://github.com/RezaNajian
#
# ==============================================================================

# ------------------------------------------------------------------------------
# Imports
# ------------------------------------------------------------------------------
from KratosMultiphysics import _ImportApplication
from KratosOptimizationApplication import *
application = KratosOptimizationApplication()
application_name = "KratosOptimizationApplication"

_ImportApplication(application, application_name)