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
class ResponseLogger():

    # --------------------------------------------------------------------------
    def InitializeLogging( self ):
        raise RuntimeError("Response logger base class is called. Please check your implementation of the function >> InitializeLogging << .")

    # --------------------------------------------------------------------------
    def LogCurrentResponses( self, optimizationIteration ):
        raise RuntimeError("Response logger base class is called. Please check your implementation of the function >> LogCurrentResponses << .")

    # --------------------------------------------------------------------------
    def FinalizeLogging( self ):
        raise RuntimeError("Response logger base class is called. Please check your implementation of the function >> FinalizeLogging << .")

    # --------------------------------------------------------------------------
    def GetValue( self, variableKey ):
        raise RuntimeError("Response logger base class is called. Please check your implementation of the function >> GetValue << .")

# ==============================================================================
