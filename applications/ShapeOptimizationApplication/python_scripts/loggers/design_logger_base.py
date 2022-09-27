# ==============================================================================
#  KratosShapeOptimizationApplication
#
#  License:         BSD License
#                   license: ShapeOptimizationApplication/license.txt
#
#  Main authors:    Baumgaertner Daniel, https://github.com/dbaumgaertner
#                   Geiser Armin, https://github.com/armingeiser
#
# ==============================================================================


# ==============================================================================
class DesignLogger:

    # --------------------------------------------------------------------------
    def InitializeLogging( self ):
        raise RuntimeError("Design logger base class is called. Please check your implementation of the function >> InitializeLogging << .")

    # --------------------------------------------------------------------------
    def LogCurrentDesign( self, optimizationIteration ):
        raise RuntimeError("Design logger base class is called. Please check your implementation of the function >> LogCurrentDesign << .")

    #---------------------------------------------------------------------------
    def FinalizeLogging( self ):
        raise RuntimeError("Design logger base class is called. Please check your implementation of the function >> FinalizeLogging << .")

# ==============================================================================
