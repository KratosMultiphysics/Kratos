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

# Import logger classes
from design_logger_gid import DesignLoggerGID
from response_logger_steepest_descent import ResponseLoggerSteepestDescent
from response_logger_penalized_projection import ResponseLoggerPenalizedProjection


# ==============================================================================
def CreateDataLogger( designSurface, communicator, timer, optimizationSettings ):
    responseLogger = createResponseLogger( communicator, timer, optimizationSettings )
    designLogger = createDesignLogger( designSurface, optimizationSettings )
    return optimizationDataLogger( responseLogger, designLogger, optimizationSettings )

# -----------------------------------------------------------------------------
def createResponseLogger( communicator, timer, optimizationSettings ):

    responseLogger = None
    optimizationAlgorithm = optimizationSettings["optimization_algorithm"]["name"].GetString()

    if optimizationAlgorithm == "steepest_descent":
        return ResponseLoggerSteepestDescent( communicator, timer, optimizationSettings )
    if optimizationAlgorithm == "penalized_projection":
        return ResponseLoggerPenalizedProjection( communicator, timer, optimizationSettings )        
    else:
        raise NameError("The following optimization algorithm not supported by the response logger (name may be a misspelling): " + optimizationAlgorithm)

# -----------------------------------------------------------------------------
def createDesignLogger( designSurface, optimizationSettings):
    designLogger = None
    outputFormatName = optimizationSettings["output"]["output_format"]["name"].GetString()
    if outputFormatName == "gid":
        return DesignLoggerGID( designSurface, optimizationSettings )
    else:
        raise NameError("The following output format is not supported by the design logger (name may be misspelled): " + outputFormatName)

# ==============================================================================
class optimizationDataLogger():

    # --------------------------------------------------------------------------
    def __init__( self, responseLogger, designLogger, optimizationSettings ):
        self.responseLogger = responseLogger
        self.designLogger = designLogger
        self.optimizationSettings = optimizationSettings
        self.createFolderToStoreOptimizationResults()     
        self.outputInformationAboutResponseFunctions()   

    # --------------------------------------------------------------------------
    def createFolderToStoreOptimizationResults ( self ):
        resultsDirectory = self.optimizationSettings["output"]["output_directory"].GetString()
        os.system( "rm -rf " + resultsDirectory )
        os.system( "mkdir -p " + resultsDirectory )          

    # --------------------------------------------------------------------------
    def outputInformationAboutResponseFunctions( self ):

        numberOfObjectives = self.optimizationSettings["objectives"].size()
        numberOfConstraints = self.optimizationSettings["constraints"].size()

        print("\n> The following objectives are defined:\n")
        for objectiveNumber in range(numberOfObjectives):
            print(self.optimizationSettings["objectives"][objectiveNumber])

        if numberOfConstraints != 0:
            print("> The following constraints are defined:\n")
            for constraintNumber in range(numberOfConstraints):
                print(self.optimizationSettings["constraints"][constraintNumber],"\n")
        else:
            print("> No constraints defined.\n")              

    # --------------------------------------------------------------------------
    def initializeDataLogging( self ):
        self.designLogger.initializeLogging()  
        self.responseLogger.initializeLogging()

    # --------------------------------------------------------------------------
    def logCurrentData( self, optimizationIteration ):    
        self.designLogger.logCurrentDesign( optimizationIteration )   
        self.responseLogger.logCurrentResponses( optimizationIteration )

    # --------------------------------------------------------------------------
    def getRelativeChangeOfObjectiveValue( self, optimizationIteration ):
        return self.responseLogger.getRelativeChangeOfObjectiveValue( optimizationIteration )

    # --------------------------------------------------------------------------
    def finalizeDataLogging( self ):
        self.designLogger.finalizeLogging()  
        self.responseLogger.finalizeLogging()

# ==============================================================================
