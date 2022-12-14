# ==============================================================================
#  KratosShapeOptimizationApplication
#
#  License:         BSD License
#                   license: ShapeOptimizationApplication/license.txt
#
#  Main authors:    Baumgaertner Daniel, https://github.com/dbaumgaertner
#
# ==============================================================================


# Additional imports
from .analyzer_base import AnalyzerBaseClass

# ==============================================================================
class EmptyAnalyzer( AnalyzerBaseClass ):
    # --------------------------------------------------------------------------
    def AnalyzeDesignAndReportToCommunicator( self, currentDesign, optimizationIteration, communicator ):
        pass

# ==============================================================================
