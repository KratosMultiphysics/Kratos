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

# importing the Kratos Library
from KratosMultiphysics import *
from KratosMultiphysics.ShapeOptimizationApplication import *

# check that KratosMultiphysics was imported in the main script
CheckForPreviousImport()

import shutil
import os

from design_logger_gid import DesignLoggerGID
from design_logger_unv import DesignLoggerUNV
from design_logger_vtk import DesignLoggerVTK

import timer_factory as timer_factory

from response_logger_steepest_descent import ResponseLoggerSteepestDescent
from response_logger_penalized_projection import ResponseLoggerPenalizedProjection

# ==============================================================================
def CreateDataLogger( OptimizationModelPart, DesignSurface, Communicator, OptimizationSettings ):
    return OptimizationDataLogger( OptimizationModelPart, DesignSurface, Communicator, OptimizationSettings )

# ==============================================================================
class OptimizationDataLogger():

    # --------------------------------------------------------------------------
    def __init__( self, OptimizationModelPart, DesignSurface, Communicator, OptimizationSettings ):
        self.OptimizationModelPart = OptimizationModelPart
        self.DesignSurface = DesignSurface
        self.Communicator = Communicator
        self.OptimizationSettings = OptimizationSettings

        self.Timer = timer_factory.CreateTimer()
        self.ResponseLogger = self.__CreateResponseLogger()
        self.DesignLogger = self.__CreateDesignLogger()

        self.__CreateFolderToStoreOptimizationResults()     
        self.__OutputInformationAboutResponseFunctions()   

    # -----------------------------------------------------------------------------
    def __CreateResponseLogger( self ):
        AlgorithmName = self.OptimizationSettings["optimization_algorithm"]["name"].GetString()
        if AlgorithmName == "steepest_descent":
            return ResponseLoggerSteepestDescent( self.Communicator, self.OptimizationSettings, self.Timer )
        elif AlgorithmName == "penalized_projection":
            return ResponseLoggerPenalizedProjection( self.Communicator, self.OptimizationSettings, self.Timer )   
        else:
            raise NameError("The following optimization algorithm not supported by the response logger (name may be a misspelling): " + AlgorithmName)

    # -----------------------------------------------------------------------------
    def __CreateDesignLogger( self ):
        outputFormatName = self.OptimizationSettings["output"]["output_format"]["name"].GetString()
        if outputFormatName == "gid":
            return DesignLoggerGID( self.OptimizationModelPart, self.DesignSurface, self.OptimizationSettings )
        if outputFormatName == "unv":
            return DesignLoggerUNV( self.OptimizationModelPart, self.DesignSurface, self.OptimizationSettings )  
        if outputFormatName == "vtk":
            return DesignLoggerVTK( self.OptimizationModelPart, self.DesignSurface, self.OptimizationSettings )                
        else:
            raise NameError("The following output format is not supported by the design logger (name may be misspelled): " + outputFormatName)

    # --------------------------------------------------------------------------
    def __CreateFolderToStoreOptimizationResults ( self ):
        resultsDirectory = self.OptimizationSettings["output"]["output_directory"].GetString()
        if os.path.exists(resultsDirectory):
            shutil.rmtree(resultsDirectory)
        os.makedirs(resultsDirectory)

    # --------------------------------------------------------------------------
    def __OutputInformationAboutResponseFunctions( self ):
        numberOfObjectives = self.OptimizationSettings["objectives"].size()
        numberOfConstraints = self.OptimizationSettings["constraints"].size()

        print("\n> The following objectives are defined:\n")
        for objectiveNumber in range(numberOfObjectives):
            print(self.OptimizationSettings["objectives"][objectiveNumber])

        if numberOfConstraints != 0:
            print("> The following constraints are defined:\n")
            for constraintNumber in range(numberOfConstraints):
                print(self.OptimizationSettings["constraints"][constraintNumber],"\n")
        else:
            print("> No constraints defined.\n")              

    # --------------------------------------------------------------------------
    def InitializeDataLogging( self ):
        self.DesignLogger.InitializeLogging()  
        self.ResponseLogger.InitializeLogging()

    # --------------------------------------------------------------------------
    def LogCurrentData( self, optimizationIteration ):        
        self.DesignLogger.LogCurrentDesign( optimizationIteration )   
        self.ResponseLogger.LogCurrentResponses( optimizationIteration )

    # --------------------------------------------------------------------------
    def FinalizeDataLogging( self ):
        self.DesignLogger.FinalizeLogging()  
        self.ResponseLogger.FinalizeLogging()

    # --------------------------------------------------------------------------
    def GetValue( self, variableKey ):
        return self.ResponseLogger.GetValue( variableKey )

    # --------------------------------------------------------------------------
    def StartTimer( self ):
        return self.Timer.StartTimer()

    # --------------------------------------------------------------------------
    def GetTimeStamp( self ):
        return self.Timer.GetTimeStamp()

    # --------------------------------------------------------------------------
    def GetLapTime( self ):
        lap_time = self.Timer.GetLapTime()
        self.Timer.StartNewLap()
        return lap_time   

    # --------------------------------------------------------------------------
    def GetTotalTime( self ):
        return self.Timer.GetTotalTime()

# ==============================================================================
