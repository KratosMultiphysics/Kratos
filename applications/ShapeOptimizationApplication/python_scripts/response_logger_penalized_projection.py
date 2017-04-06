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
class ResponseLoggerPenalizedProjection( ResponseLogger ):

    # --------------------------------------------------------------------------
    def __init__( self, communicator , optimizationSettings, timer, additionalVariablesToBeLogged ):
        self.communicator = communicator
        self.optimizationSettings = optimizationSettings
        self.timer = timer
        self.stepSize = additionalVariablesToBeLogged["stepSize"]
        self.correctionScaling = additionalVariablesToBeLogged["correctionScaling"]

        self.completeResponseLogFileName = self.createCompleteResponseLogFilename( optimizationSettings )
        self.onlyObjective = self.optimizationSettings["objectives"][0]["identifier"].GetString()   
        self.onlyConstraint = self.optimizationSettings["constraints"][0]["identifier"].GetString()  
        self.typOfOnlyConstraint = self.optimizationSettings["constraints"][0]["type"].GetString()    

        self.objectiveValueHistory = {}
        self.constraintValueHistory = {}
        self.constraintReferenceValueHistory = {}
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
            row.append("\tc["+self.onlyConstraint+"]: "+self.typOfOnlyConstraint+"\t")    
            row.append("\tc["+self.onlyConstraint+"] / reference_value[%]"+"\t")        
            row.append("\tcorrection_scaling[-]\t")
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

        self.addCurrentResponseValuesToHistory()
        if self.isFirstLog():
            self.initializeChangeOfObjectiveValueHistory()
        else:
            self.addChangeOfObjectiveValueToHistory()
            
        self.printInfoAboutResponseFunction()
        self.writeDataToLogFile()

        self.previousOptimizationIteration = optimizationIteration

    # -------------------------------------------------------------------------
    def isFirstLog( self ):
        if len(self.objectiveValueHistory) < 2:
            return True
        else:
            return False

    # --------------------------------------------------------------------------
    def addCurrentResponseValuesToHistory( self ):
        self.addObjectiveValuesToHistory()
        self.addConstraintValuesToHistory()

    # --------------------------------------------------------------------------
    def addObjectiveValuesToHistory( self ):
        objectiveValue = self.communicator.getReportedFunctionValueOf ( self.onlyObjective )
        self.objectiveValueHistory[self.currentOptimizationIteration] = objectiveValue

    # --------------------------------------------------------------------------
    def addConstraintValuesToHistory( self ):
        constraintValue = self.communicator.getReportedFunctionValueOf ( self.onlyConstraint )
        constraintReferenceValue = self.communicator.getReportedFunctionReferenceValueOf ( self.onlyConstraint )
        self.constraintValueHistory[self.currentOptimizationIteration] = constraintValue
        self.constraintReferenceValueHistory[self.currentOptimizationIteration] = constraintReferenceValue

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
    def printInfoAboutResponseFunction( self ):
        objectiveValue = self.objectiveValueHistory[self.currentOptimizationIteration]
        constraintValue = self.constraintValueHistory[self.currentOptimizationIteration]
        absoluteChangeOfObjectiveValue = self.absoluteChangeOfObjectiveValueHistory[self.currentOptimizationIteration]
        relativeChangeOfObjectiveValue = self.relativeChangeOfObjectiveValueHistory[self.currentOptimizationIteration]        
        print("\n> Current value of objective function = ", round(objectiveValue,12))
        print("> Absolut change of objective function = ",round(absoluteChangeOfObjectiveValue,4)," [%]")
        print("> Relative change of objective function = ",round(relativeChangeOfObjectiveValue,4)," [%]")         
        print("\n> Current value of constraint function = ", round(constraintValue,12))

    # --------------------------------------------------------------------------
    def writeDataToLogFile( self ):

        objectiveValue = self.objectiveValueHistory[self.currentOptimizationIteration]
        constraintValue = self.constraintValueHistory[self.currentOptimizationIteration]
        constraintReferenceValue = self.constraintReferenceValueHistory[self.currentOptimizationIteration]
        absoluteChangeOfObjectiveValue = self.absoluteChangeOfObjectiveValueHistory[self.currentOptimizationIteration]
        relativeChangeOfObjectiveValue = self.relativeChangeOfObjectiveValueHistory[self.currentOptimizationIteration]             

        with open(self.completeResponseLogFileName, 'a') as csvfile:
            historyWriter = csv.writer(csvfile, delimiter=',',quotechar='|',quoting=csv.QUOTE_MINIMAL)
            row = []
            row.append(str(self.currentOptimizationIteration)+"\t")
            row.append("\t"+str("%.12f"%(objectiveValue))+"\t")
            row.append("\t"+str("%.2f"%(absoluteChangeOfObjectiveValue))+"\t")
            row.append("\t"+str("%.6f"%(relativeChangeOfObjectiveValue))+"\t")
            row.append("\t"+str("%.12f"%(constraintValue))+"\t")
            if not constraintReferenceValue:
                row.append("\t"+str("-\t"))
            else: 
                percentageOfReference = 100*(constraintValue / constraintReferenceValue)
                row.append("\t"+str("%.6f"%(percentageOfReference)))
            row.append("\t"+str("%.6f"%(self.correctionScaling))+"\t")
            row.append("\t"+str(self.stepSize)+"\t")
            row.append("\t"+str("%.1f"%(self.timer.getLapTime()))+"\t")
            row.append("\t"+str("%.1f"%(self.timer.getTotalTime()))+"\t")
            row.append("\t"+str(self.timer.getTimeStamp()))
            historyWriter.writerow(row)   

    # --------------------------------------------------------------------------
    def getValue( self, variableKey ):
        if variableKey=="RELATIVE_CHANGE_OF_OBJECTIVE_VALUE":
            if self.isFirstLog():
                raise RuntimeError("Relative change of objective function can not be computed since only one logged value is existing!")
            else:
                return self.relativeChangeOfObjectiveValueHistory[self.currentOptimizationIteration]
        else:
            raise NameError("Value with the following variable key not defined in response_logger_penalized_projection.py: " + variableKey)

    # --------------------------------------------------------------------------
    def finalizeLogging( self ):      
        pass # No finalization necessary here

# ==============================================================================
