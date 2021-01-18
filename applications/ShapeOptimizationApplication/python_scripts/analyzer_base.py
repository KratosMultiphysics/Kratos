# ==============================================================================
#  KratosShapeOptimizationApplication
#
#  License:         BSD License
#                   license: ShapeOptimizationApplication/license.txt
#
#  Main authors:    Baumgaertner Daniel, https://github.com/dbaumgaertner
#
# ==============================================================================

# Making KratosMultiphysics backward compatible with python 2.6 and 2.7
from __future__ import print_function, absolute_import, division

# ==============================================================================
class AnalyzerBaseClass:
    # --------------------------------------------------------------------------
    def __init__(self):
        pass

    # --------------------------------------------------------------------------
    def InitializeBeforeOptimizationLoop(self):
        pass

    # --------------------------------------------------------------------------
    def AnalyzeDesignAndReportToCommunicator(self, current_design, optimization_iteration, communicator):
        raise RuntimeError("Analyzer base class is called. Please check your implementation of the function >> AnalyzeDesignAndReportToCommunicator << .")

    # --------------------------------------------------------------------------
    def FinalizeAfterOptimizationLoop(self):
        pass

# ==============================================================================
