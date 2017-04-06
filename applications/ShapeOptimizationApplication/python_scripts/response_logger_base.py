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
class ResponseLogger():

    # --------------------------------------------------------------------------
    def __init__( self, communicator, timer, optimizationSettings ):
        raise RuntimeError("Response logger base class is called. Please check your implementation of the function >> __init__ << .")

    # --------------------------------------------------------------------------
    def initializeLogging( self ):
        raise RuntimeError("Response logger base class is called. Please check your implementation of the function >> initializeLogging << .")

    # --------------------------------------------------------------------------
    def logCurrentResponses( self, optimizationIteration ):
        raise RuntimeError("Response logger base class is called. Please check your implementation of the function >> logCurrentResponses << .")

    # --------------------------------------------------------------------------
    def getRelativeChangeOfObjectiveValue( self, optimizationIteration ):
        raise RuntimeError("Response logger base class is called. Please check your implementation of the function >> getRelativeChangeOfObjectiveValue << .")

    # --------------------------------------------------------------------------
    def finalizeLogging( self ):      
        raise RuntimeError("Response logger base class is called. Please check your implementation of the function >> finalizeLogging << .")

# ==============================================================================
