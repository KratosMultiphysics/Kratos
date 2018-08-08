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

# Making KratosMultiphysics backward compatible with python 2.6 and 2.7
from __future__ import print_function, absolute_import, division

# ==============================================================================
class ValueLogger():

    # --------------------------------------------------------------------------
    def InitializeLogging( self ):
        raise RuntimeError("Value logger base class is called. Please check your implementation of the function >> InitializeLogging << .")

    # --------------------------------------------------------------------------
    def LogCurrentValues( self, optimizationIteration ):
        raise RuntimeError("Value logger base class is called. Please check your implementation of the function >> LogCurrentValues << .")

    # --------------------------------------------------------------------------
    def FinalizeLogging( self ):
        raise RuntimeError("Value logger base class is called. Please check your implementation of the function >> FinalizeLogging << .")

    # --------------------------------------------------------------------------
    def GetHistoryOfValues( self ):
        raise RuntimeError("Value logger base class is called. Please check your implementation of the function >> GetHistoryOfValues << .")


# ==============================================================================
