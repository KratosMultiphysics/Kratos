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
from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
from KratosMultiphysics import _ImportApplication
from KratosTopologyOptimizationApplication import *
application = KratosTopologyOptimizationApplication()
application_name = "KratosTopologyOptimizationApplication"

_ImportApplication(application, application_name)
