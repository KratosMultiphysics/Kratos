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

# Import logger base classes
from response_logger_base import ResponseLogger

# Import additional libraries
import csv

# ==============================================================================
class ResponseLoggerPenalizedProjection( ResponseLogger ):

    # --------------------------------------------------------------------------
    def __init__( self, communicator, optimizationSettings, timer ):
        self.communicator = communicator
        self.optimizationSettings = optimizationSettings
        self.timer = timer

        self.onlyObjective = self.optimizationSettings["objectives"][0]["identifier"].GetString()
        self.onlyConstraint = self.optimizationSettings["constraints"][0]["identifier"].GetString()
        self.typOfOnlyConstraint = self.optimizationSettings["constraints"][0]["type"].GetString()

        self.completeResponseLogFileName = self.__CreateCompleteResponseLogFilename( optimizationSettings )

        self.objectiveValueHistory = {}
        self.constraintValueHistory = {}
        self.objectiveReferenceValue = None
        self.absoluteChangeOfObjectiveValueHistory = {}
        self.relativeChangeOfObjectiveValueHistory = {}

        self.currentIteration = 0
        self.previousIteration = 0
        self.initialIteration = 0

    # --------------------------------------------------------------------------
    def __CreateCompleteResponseLogFilename( self, optimizationSettings ):
        resultsDirectory = optimizationSettings["output"]["output_directory"].GetString()
        responseLogFilename = optimizationSettings["output"]["response_log_filename"].GetString()
        completeResponseLogFilename = resultsDirectory+"/"+responseLogFilename+".csv"
        return completeResponseLogFilename

    # --------------------------------------------------------------------------
    def InitializeLogging( self ):
        with open(self.completeResponseLogFileName, 'w') as csvfile:
            historyWriter = csv.writer(csvfile, delimiter=',',quotechar='|',quoting=csv.QUOTE_MINIMAL)
            row = []
            row.append("{:<6s}".format("itr"))
            row.append("{:>20s}".format("f"))
            row.append("{:>12s}".format("df_abs[%]"))
            row.append("{:>12s}".format("df_rel[%]"))
            row.append("{:>20s}".format("c["+self.onlyConstraint+"]: "+self.typOfOnlyConstraint))
            row.append("{:>20s}".format("c["+self.onlyConstraint+"]_ref"))
            row.append("{:>13s}".format("c_scaling[-]"))
            row.append("{:>13s}".format("step_size[-]"))
            row.append("{:>12s}".format("t_itr[s]"))
            row.append("{:>16s}".format("t_total[s]"))
            row.append("{:>25s}".format("time_stamp"))
            historyWriter.writerow(row)

    # --------------------------------------------------------------------------
    def LogCurrentResponses( self, optimizationIteration ):

        self.currentIteration = optimizationIteration

        if self.__IsFirstLog():
            self.initialIteration = optimizationIteration
            self.__AddResponseValuesToHistory()
            self.__InitializeChangeOfResponseValuesHistory()
        else:
            self.__AddResponseValuesToHistory()
            self.__AddChangeOfResponseValuesToHistory()
        self.__PrintInfoAboutResponseFunctionValues()
        self.__WriteDataToLogFile()

        self.previousIteration = optimizationIteration

    # -------------------------------------------------------------------------
    def __IsFirstLog( self ):
        if len(self.objectiveValueHistory) == 0:
            return True
        else:
            return False

    # --------------------------------------------------------------------------
    def __AddResponseValuesToHistory( self ):
        objectiveValue = self.communicator.getFunctionValue( self.onlyObjective )
        self.objectiveValueHistory[self.currentIteration] = objectiveValue

        constraintValue = self.communicator.getFunctionValue( self.onlyConstraint )
        self.constraintValueHistory[self.currentIteration] = constraintValue

    # --------------------------------------------------------------------------
    def __InitializeChangeOfResponseValuesHistory( self ):
        self.objectiveReferenceValue = self.objectiveValueHistory[self.initialIteration]
        if abs(self.objectiveReferenceValue)<1e-12:
            print("\n> WARNING: Objective reference value < 1e-12!!:")
            print("> WARNING: Standard reference value of 1 is assumed.")
            self.objectiveReferenceValue = 1.0

        self.absoluteChangeOfObjectiveValueHistory[self.currentIteration] = 0.0
        self.relativeChangeOfObjectiveValueHistory[self.currentIteration] = 0.0

    # --------------------------------------------------------------------------
    def __AddChangeOfResponseValuesToHistory( self ):
        objectiveValue = self.objectiveValueHistory[self.currentIteration]
        previousObjectiveValue = self.objectiveValueHistory[self.previousIteration]
        initialObjectiveValue = self.objectiveValueHistory[self.initialIteration]

        constraintValue = self.constraintValueHistory[self.currentIteration]
        previousConstraintValue = self.constraintValueHistory[self.previousIteration]
        initialConstraintValue = self.constraintValueHistory[self.initialIteration]

        self.absoluteChangeOfObjectiveValueHistory[self.currentIteration] = 100*(objectiveValue-initialObjectiveValue) / abs(self.objectiveReferenceValue)
        self.relativeChangeOfObjectiveValueHistory[self.currentIteration] = 100*(objectiveValue-previousObjectiveValue) / abs(self.objectiveReferenceValue)

    # --------------------------------------------------------------------------
    def __PrintInfoAboutResponseFunctionValues( self ):
        objectiveValue = self.objectiveValueHistory[self.currentIteration]
        constraintValue = self.constraintValueHistory[self.currentIteration]

        absoluteChangeOfObjectiveValue = self.absoluteChangeOfObjectiveValueHistory[self.currentIteration]
        relativeChangeOfObjectiveValue = self.relativeChangeOfObjectiveValueHistory[self.currentIteration]

        print("\n> Current value of objective function = ", round(objectiveValue,12))
        print("> Absolut change of objective function = ",round(absoluteChangeOfObjectiveValue,4)," [%]")
        print("> Relative change of objective function = ",round(relativeChangeOfObjectiveValue,4)," [%]")
        print("\n> Current value of constraint function = ", round(constraintValue,12))

    # --------------------------------------------------------------------------
    def __WriteDataToLogFile( self ):
        objectiveValue = self.objectiveValueHistory[self.currentIteration]
        constraintValue = self.constraintValueHistory[self.currentIteration]
        constraintReferenceValue = self.communicator.getReferenceValue( self.onlyConstraint )

        absoluteChangeOfObjectiveValue = self.absoluteChangeOfObjectiveValueHistory[self.currentIteration]
        relativeChangeOfObjectiveValue = self.relativeChangeOfObjectiveValueHistory[self.currentIteration]

        with open(self.completeResponseLogFileName, 'a') as csvfile:
            historyWriter = csv.writer(csvfile, delimiter=',',quotechar='|',quoting=csv.QUOTE_MINIMAL)
            row = []
            row.append("{:<6s}".format(str(self.currentIteration)))
            row.append(str("{:>20f}".format(objectiveValue)))
            row.append(str("{:>12f}".format(absoluteChangeOfObjectiveValue)))
            row.append(str("{:>12f}".format(relativeChangeOfObjectiveValue)))
            row.append(str("{:>20f}".format(constraintValue)))
            row.append(str("{:>20f}".format(constraintReferenceValue)))
            row.append(str("{:>13f}".format(self.optimizationSettings["optimization_algorithm"]["correction_scaling"].GetDouble())))
            row.append(str("{:>13f}".format(self.optimizationSettings["line_search"]["step_size"].GetDouble())))
            row.append(str("{:>12f}".format(self.timer.GetLapTime())))
            row.append(str("{:>16f}".format(self.timer.GetTotalTime())))
            row.append("{:>25}".format(self.timer.GetTimeStamp()))
            historyWriter.writerow(row)

    # --------------------------------------------------------------------------
    def FinalizeLogging( self ):
        pass # No finalization necessary here

    # --------------------------------------------------------------------------
    def GetValue( self, variableKey ):
        if variableKey=="RELATIVE_CHANGE_OF_OBJECTIVE_VALUE":
            if self.__IsFirstLog():
                raise RuntimeError("Relative change of objective function can not be computed since only one logged value is existing!")
            else:
                return self.relativeChangeOfObjectiveValueHistory[self.currentIteration]
        else:
            raise NameError("Value with the following variable key not defined in response_logger_penalized_projection.py: " + variableKey)

# ==============================================================================
