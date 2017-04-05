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
from optimization_data_logger.response_logger_base import ResponseLogger

# Import additional libraries
import csv
import time

# ==============================================================================
class ResponseLoggerSteepestDescent():

    # --------------------------------------------------------------------------
    def __init__( self, optimizationSettings ):
        self.optimizationSettings = optimizationSettings
        self.initialObjectiveValue = None
        self.previousObjectiveValue = None
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
    def logCurrentResponses( self, optimizationIteration, communicator ):
        
        # Get Id of objective (right now only one objective is considered!!)
        onlyObjectiveFunction = self.optimizationSettings["objectives"][0]["identifier"].GetString()
        valueOfObjectiveFunction = communicator.getReportedFunctionValueOf ( onlyObjectiveFunction )
        print("\n> Current value of objective function = ",valueOfObjectiveFunction)
        if optimizationIteration == 1:
            self.saveInitialValues( valueOfObjectiveFunction )   

        absoluteChangeOfObjectiveValue = 0.0
        relativeChangeOfObjectiveValue = 0.0
        if optimizationIteration > 1:
            absoluteChangeOfObjectiveValue = 100*(valueOfObjectiveFunction-self.initialValueOfObjectiveFunction) / self.initialValueOfObjectiveFunction
            relativeChangeOfObjectiveValue = 100*(valueOfObjectiveFunction-self.previousValueOfObjectiveFunction) / self.initialValueOfObjectiveFunction
            print("> Absolut change of objective function = ",round(absoluteChangeOfObjectiveValue,6)," [%]")
            print("> Relative change of objective function = ",round(relativeChangeOfObjectiveValue,6)," [%]") 


        runTimeOptimizationStep = 666666666666666666666# round(time.time() - timeAtStartOfCurrentOptimizationStep,2)
        runTimeOptimization = 6666666666666666666666666#round(time.time() - self.timeAtStartOfOptimization,2)
        with open(self.completeResponseLogFileName, 'a') as csvfile:
            historyWriter = csv.writer(csvfile, delimiter=',',quotechar='|',quoting=csv.QUOTE_MINIMAL)
            row = []
            row.append(str(optimizationIteration)+"\t")
            row.append("\t"+str("%.12f"%(valueOfObjectiveFunction))+"\t")
            row.append("\t"+str("%.2f"%(absoluteChangeOfObjectiveValue))+"\t")
            row.append("\t"+str("%.6f"%(relativeChangeOfObjectiveValue))+"\t")
            row.append("\t"+str(self.optimizationSettings["line_search"]["step_size"].GetDouble())+"\t")
            row.append("\t"+str("%.1f"%(runTimeOptimizationStep))+"\t")
            row.append("\t"+str("%.1f"%(runTimeOptimization))+"\t")
            row.append("\t"+str(time.ctime()))
            historyWriter.writerow(row)  
            
        self.saveValuesForNextOptimizationIteration( valueOfObjectiveFunction ) 

    # --------------------------------------------------------------------------
    def saveInitialValues( self, valueOfObjectiveFunction ):
        self.initialValueOfObjectiveFunction = valueOfObjectiveFunction

    # --------------------------------------------------------------------------
    def saveValuesForNextOptimizationIteration( self, valueOfObjectiveFunction ):
        self.previousValueOfObjectiveFunction = valueOfObjectiveFunction        

    # --------------------------------------------------------------------------
    def finalizeLogging( self ):      
        pass # No finalization necessary here

# ==============================================================================
