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
    def __init__( self, communicator, timer, optimizationSettings ):
        self.communicator = communicator
        self.timer = timer
        self.optimizationSettings = optimizationSettings
        self.objectiveValueHistory = {}
        self.currentOptimizationIteration = 0
        self.previousOptimizationIteration = 0
        self.initialOptimizationIteration = 0
        self.absoluteChangeOfObjectiveValue = 0.0
        self.relativeChangeOfObjectiveValue = 0.0
        self.completeResponseLogFileName = self.createCompleteResponseLogFilename( optimizationSettings )

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
        if self.isFirstLog():
            self.initialOptimizationIteration = optimizationIteration
            self.addCurrentResponseFunctionValuesToHistory()
        else:
            self.addCurrentResponseFunctionValuesToHistory()
            self.evaluateChangeOfResponseFunctionValues()
        self.writeDataToLogFile()
        self.previousOptimizationIteration = optimizationIteration

    # -------------------------------------------------------------------------
    def isFirstLog( self ):
        if len(self.objectiveValueHistory) == 0:
            return True
        else:
            return False

    # --------------------------------------------------------------------------
    def addCurrentResponseFunctionValuesToHistory( self ):
        onlyObjective = self.optimizationSettings["objectives"][0]["identifier"].GetString()
        objectiveValue = self.communicator.getReportedFunctionValueOf ( onlyObjective )
        self.objectiveValueHistory[self.currentOptimizationIteration] = objectiveValue
        print("\n> Current value of objective function = ", round(objectiveValue,12))

    # --------------------------------------------------------------------------
    def evaluateChangeOfResponseFunctionValues( self ):
        objectiveValue = self.objectiveValueHistory[self.currentOptimizationIteration]
        previousObjectiveValue = self.objectiveValueHistory[self.previousOptimizationIteration]
        initialObjectiveValue = self.objectiveValueHistory[self.initialOptimizationIteration]
        self.absoluteChangeOfObjectiveValue = 100*(objectiveValue-initialObjectiveValue) / initialObjectiveValue
        self.relativeChangeOfObjectiveValue = 100*(objectiveValue-previousObjectiveValue) / initialObjectiveValue
        print("> Absolut change of objective function = ",round(self.absoluteChangeOfObjectiveValue,4)," [%]")
        print("> Relative change of objective function = ",round(self.relativeChangeOfObjectiveValue,4)," [%]") 

    # --------------------------------------------------------------------------
    def writeDataToLogFile( self ):
        objectiveValue = self.objectiveValueHistory[self.currentOptimizationIteration]
        with open(self.completeResponseLogFileName, 'a') as csvfile:
            historyWriter = csv.writer(csvfile, delimiter=',',quotechar='|',quoting=csv.QUOTE_MINIMAL)
            row = []
            row.append(str(self.currentOptimizationIteration)+"\t")
            row.append("\t"+str("%.12f"%(objectiveValue))+"\t")
            row.append("\t"+str("%.2f"%(self.absoluteChangeOfObjectiveValue))+"\t")
            row.append("\t"+str("%.6f"%(self.relativeChangeOfObjectiveValue))+"\t")
            row.append("\t"+str(self.optimizationSettings["line_search"]["step_size"].GetDouble())+"\t")
            row.append("\t"+str("%.1f"%(self.timer.getLapTime()))+"\t")
            row.append("\t"+str("%.1f"%(self.timer.getTotalTime()))+"\t")
            row.append("\t"+str(self.timer.getTimeStamp()))
            historyWriter.writerow(row)       

    # --------------------------------------------------------------------------
    def getRelativeChangeOfObjectiveValue( self ):
        if self.isFirstLog():
            raise RuntimeError("Relative change of objective function can not be computed since only one logged value is existing!")
        else:
            return self.relativeChangeOfObjectiveValue

    # --------------------------------------------------------------------------
    def finalizeLogging( self ):      
        pass # No finalization necessary here

# ==============================================================================
