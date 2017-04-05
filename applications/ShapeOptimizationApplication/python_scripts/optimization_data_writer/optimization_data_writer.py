# ==============================================================================
#  KratosShapeOptimizationApplication
#
#  License:         BSD License
#                   license: ShapeOptimizationApplication/license.txt
#
#  Main authors:    Baumgärtner Daniel, https://github.com/dbaumgaertner
#                   Geiser Armin, https://github.com/armingeiser
#
# ==============================================================================

# Making KratosMultiphysics backward compatible with python 2.6 and 2.7
from __future__ import print_function, absolute_import, division

# importing the Kratos Library
from KratosMultiphysics import *
from KratosMultiphysics.ShapeOptimizationApplication import *

# check that KratosMultiphysics was imported in the main script
CheckForPreviousImport()

# ==============================================================================
def optimizationDataWriter():

    # --------------------------------------------------------------------------
    def __init__( self ):

    # --------------------------------------------------------------------------
    def initializeDesignOutput():
        raise RuntimeError("opimizationDataWriter base class is called. Please check your implementation of the function >> initializeDesignOutput << .")

    # --------------------------------------------------------------------------
    def initializeOptimizationLog():
        raise RuntimeError("opimizationDataWriter base class is called. Please check your implementation of the function >> initializeOptimizationLog << .")

    # --------------------------------------------------------------------------
    def writeCurrentDesign():
        raise RuntimeError("opimizationDataWriter base class is called. Please check your implementation of the function >> writeCurrentDesign << .")

    # --------------------------------------------------------------------------
    def logCurrentOptimizationStep():
        raise RuntimeError("opimizationDataWriter base class is called. Please check your implementation of the function >> logCurrentOptimizationStep << .")

    # --------------------------------------------------------------------------
    def finalizeDesignOutput():
        raise RuntimeError("opimizationDataWriter base class is called. Please check your implementation of the function >> finalizeDesignOutput << .")

    # --------------------------------------------------------------------------
    def fianlizeOptimizationLog():    
        raise RuntimeError("opimizationDataWriter base class is called. Please check your implementation of the function >> fianlizeOptimizationLog << .")

# ==============================================================================
