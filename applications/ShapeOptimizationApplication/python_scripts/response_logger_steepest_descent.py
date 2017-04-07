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

# Import logger base classes
from response_logger_base import ResponseLogger

# Import additional libraries
import csv

# ==============================================================================
class ResponseLoggerSteepestDescent( ResponseLogger ):

    # --------------------------------------------------------------------------
    def __init__( self, communicator, optimizationSettings, timer, specificVariablesToBeLogged ):

        self.communicator = communicator
        self.optimizationSettings = optimizationSettings
        self.timer = timer
        self.specificVariablesToBeLogged = specificVariablesToBeLogged

        self.onlyObjective = self.optimizationSettings["objectives"][0]["identifier"].GetString()   

        self.completeResponseLogFileName = self.createCompleteResponseLogFilename( optimizationSettings )

        self.objectiveValueHistory = {}
        self.absoluteChangeOfObjectiveValueHistory = {}
        self.relativeChangeOfObjectiveValueHistory = {}

        self.currentOptimizationIteration = 0
        self.previousOptimizationIteration = 0
        self.initialOptimizationIteration = 0

    # --------------------------------------------------------------------------
    def createCompleteResponseLogFilename( self, optimizationSettings ):
        resultsDirectory = optimizationSettings["output"]["output_directory"].GetString()
        responseLogFilename = optimizationSettings["output"]["response_log_filename"].GetString()
        completeResponseLogFilename = resultsDirectory+"/"+responseLogFilename+".csv"
        return completeResponseLogFilename     

    # --------------------------------------------------------------------------
    def initializeLogging( self ):
        with open(self.completeResponseLogFileName, 'w') as csvfile:
            historyWriter = csv.writer(csvfile, delimiter=',',quotechar='|',quoting=csv.QUOTE_MINIMAL)
            row = []
            row.append("itr\t")
            row.append("\tf\t")
            row.append("\tdf_absolute[%]\t")
            row.append("\tdf_relative[%]\t")
            row.append("\tstep_size[-]\t")
            row.append("\tt_iteration[s]\t")
            row.append("\tt_total[s]") 
            row.append("\ttime_stamp") 
            historyWriter.writerow(row)    

    # --------------------------------------------------------------------------
    def logCurrentResponses( self, optimizationIteration ):
        self.currentOptimizationIteration = optimizationIteration

        self.addCurrentObjectiveValueToHistory()
        if self.isFirstLog():
            self.initialOptimizationIteration = optimizationIteration        
            self.initializeChangeOfObjectiveValueHistory()
        else:
            self.addChangeOfObjectiveValueToHistory()
        self.printInfoAboutResponseFunctionValues()
        self.writeDataToLogFile()
        
        self.previousOptimizationIteration = optimizationIteration

    # -------------------------------------------------------------------------
    def isFirstLog( self ):
        if len(self.objectiveValueHistory) < 2:
            return True
        else:
            return False

    # --------------------------------------------------------------------------
    def addCurrentObjectiveValueToHistory( self ):
        objectiveValue = self.communicator.getReportedFunctionValueOf ( self.onlyObjective )
        self.objectiveValueHistory[self.currentOptimizationIteration] = objectiveValue
        
    # --------------------------------------------------------------------------
    def initializeChangeOfObjectiveValueHistory( self ):
        self.absoluteChangeOfObjectiveValueHistory[self.currentOptimizationIteration] = 0.0
        self.relativeChangeOfObjectiveValueHistory[self.currentOptimizationIteration] = 0.0

    # --------------------------------------------------------------------------
    def addChangeOfObjectiveValueToHistory( self ):
        objectiveValue = self.objectiveValueHistory[self.currentOptimizationIteration]
        previousObjectiveValue = self.objectiveValueHistory[self.previousOptimizationIteration]
        initialObjectiveValue = self.objectiveValueHistory[self.initialOptimizationIteration]
        self.absoluteChangeOfObjectiveValueHistory[self.currentOptimizationIteration] = 100*(objectiveValue-initialObjectiveValue) / initialObjectiveValue
        self.relativeChangeOfObjectiveValueHistory[self.currentOptimizationIteration] = 100*(objectiveValue-previousObjectiveValue) / initialObjectiveValue

    # --------------------------------------------------------------------------
    def printInfoAboutResponseFunctionValues( self ):
        objectiveValue = self.objectiveValueHistory[self.currentOptimizationIteration]
        absoluteChangeOfObjectiveValue = self.absoluteChangeOfObjectiveValueHistory[self.currentOptimizationIteration]
        relativeChangeOfObjectiveValue = self.relativeChangeOfObjectiveValueHistory[self.currentOptimizationIteration]        
        print("\n> Current value of objective function = ", round(objectiveValue,12))
        print("> Absolut change of objective function = ",round(absoluteChangeOfObjectiveValue,4)," [%]")
        print("> Relative change of objective function = ",round(relativeChangeOfObjectiveValue,4)," [%]")         

    # --------------------------------------------------------------------------
    def writeDataToLogFile( self ):

        objectiveValue = self.objectiveValueHistory[self.currentOptimizationIteration]
        absoluteChangeOfObjectiveValue = self.absoluteChangeOfObjectiveValueHistory[self.currentOptimizationIteration]
        relativeChangeOfObjectiveValue = self.relativeChangeOfObjectiveValueHistory[self.currentOptimizationIteration]  
        stepSize = self.specificVariablesToBeLogged["stepSize"]        
        
        with open(self.completeResponseLogFileName, 'a') as csvfile:
            historyWriter = csv.writer(csvfile, delimiter=',',quotechar='|',quoting=csv.QUOTE_MINIMAL)
            row = []
            row.append(str(self.currentOptimizationIteration)+"\t")
            row.append("\t"+str("%.12f"%(objectiveValue))+"\t")
            row.append("\t"+str("%.2f"%(absoluteChangeOfObjectiveValue))+"\t")
            row.append("\t"+str("%.6f"%(relativeChangeOfObjectiveValue))+"\t")
            row.append("\t"+str(stepSize)+"\t")
            row.append("\t"+str("%.1f"%(self.timer.getLapTime()))+"\t")
            row.append("\t"+str("%.1f"%(self.timer.getTotalTime()))+"\t")
            row.append("\t"+str(self.timer.getTimeStamp()))
            historyWriter.writerow(row)       

    # --------------------------------------------------------------------------
    def finalizeLogging( self ):      
        pass # No finalization necessary here

    # --------------------------------------------------------------------------
    def getValue( self, variableKey ):
        
        if variableKey=="RELATIVE_CHANGE_OF_OBJECTIVE_VALUE":
            if self.isFirstLog():
                raise RuntimeError("Relative change of objective function can not be computed since only one logged value is existing!")
            else:
                return self.relativeChangeOfObjectiveValueHistory[self.currentOptimizationIteration]
        else:
            raise NameError("Value with the following variable key not defined in response_logger_penalized_projection.py: " + variableKey)

# ==============================================================================
