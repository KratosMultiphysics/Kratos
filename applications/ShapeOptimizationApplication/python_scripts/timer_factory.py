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

# Import additional libraries
import time

# ==============================================================================
def CreateTimer():
    return Timer()

# ==============================================================================
class Timer:

    # --------------------------------------------------------------------------
    def __init__( self ):
        self.precision = 3
        self.startTimeGlobal = None
        self.startTimeLap = None
        self.lapTime = None
    # --------------------------------------------------------------------------
    def startTimer( self ):
        self.startTimeGlobal = time.time()
        self.startTimeLap = time.time()

    # --------------------------------------------------------------------------
    def getLapTime( self ):
        lapTime = time.time() - self.startTimeLap
        return round( lapTime, self.precision )

    # --------------------------------------------------------------------------
    def resetLapTime( self ):    
        self.startTimeLap = time.time()    

    #---------------------------------------------------------------------------
    def getTotalTime( self ):
        totalTime = time.time() - self.startTimeGlobal
        return round( totalTime, self.precision )

    # --------------------------------------------------------------------------
    def getTimeStamp( self ):
        return time.ctime()

# ==============================================================================
