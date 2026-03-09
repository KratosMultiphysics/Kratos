# ==============================================================================
#  KratosOptimizationApplication
#
#  License:         BSD License
#                   license: OptimizationApplication/license.txt
#
#  Main authors:    Reza Najian Asl, https://github.com/RezaNajian
#
# ==============================================================================

# Import additional libraries
import time

# ==============================================================================
class Timer:
    # --------------------------------------------------------------------------
    def __init__( self ):
        self.Precision = 3
        self.StartGlobalTime = None
        self.StartLapTime = None
        self.LapTime = None
    # --------------------------------------------------------------------------
    def StartTimer( self ):
        self.StartGlobalTime = time.time()
        self.StartLapTime = time.time()

    # --------------------------------------------------------------------------
    def GetLapTime( self ):
        LapTime = time.time() - self.StartLapTime
        return round( LapTime, self.Precision )

    # --------------------------------------------------------------------------
    def StartNewLap( self ):
        self.StartLapTime = time.time()

    #---------------------------------------------------------------------------
    def GetTotalTime( self ):
        totalTime = time.time() - self.StartGlobalTime
        return round( totalTime, self.Precision )

    # --------------------------------------------------------------------------
    def GetTimeStamp( self ):
        return time.ctime()

# ==============================================================================
