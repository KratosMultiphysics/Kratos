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
