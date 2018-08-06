# ==============================================================================
#  KratosShapeOptimizationApplication
#
#  License:         BSD License
#                   license: ShapeOptimizationApplication/license.txt
#
#  Main authors:    Baumgaertner Daniel, https://github.com/dbaumgaertner
#
# ==============================================================================

# Making KratosMultiphysics backward compatible with python 2.6 and 2.7
from __future__ import print_function, absolute_import, division

# Kratos Core and Apps
from KratosMultiphysics import *
from KratosMultiphysics.ShapeOptimizationApplication import *

# Import logger base classes
from response_logger_base import ResponseLogger

# Import additional libraries
import csv
from custom_timer import Timer

# ==============================================================================
class ResponseLoggerTrustRegion( ResponseLogger ):

    # --------------------------------------------------------------------------
    def __init__( self, optimization_settings ):
        self.optimization_settings = optimization_settings

        self.specified_objectives = optimization_settings["objectives"]
        self.specified_constraints = optimization_settings["constraints"]

        self.complete_response_log_filename = self.__CreateCompleteResponseLogFilename( optimization_settings )

        self.objective_history = {}
        self.constraints_history = {}
        for itr in range(self.specified_constraints.size()):
            Ci_id = self.specified_constraints[itr]["identifier"].GetString()
            self.constraints_history[Ci_id] = {}

        self.objective_output_reference = None
        self.constraints_output_reference = {}

        self.abs_change_of_objective = 0.0
        self.rel_change_of_objective = 0.0

        self.current_iteration = 0
        self.previous_iteration = 0
        self.initial_iteration = 0

    # --------------------------------------------------------------------------
    def InitializeLogging( self ):
        with open(self.complete_response_log_filename, 'w') as csvfile:
            historyWriter = csv.writer(csvfile, delimiter=',',quotechar='|',quoting=csv.QUOTE_MINIMAL)
            row = []
            row.append("itr")
            row.append("\tf")
            row.append("\tdf_abs[%]")
            row.append("\tdf_rel[%]")

            for itr in range(self.specified_constraints.size()):
                con = self.specified_constraints[itr]
                con_type = con["type"].GetString()
                row.append("\tc["+str(itr)+"]: "+con_type)
                row.append("\tc["+str(itr)+"]_ref")

            row.append("\tstep_length[-]")
            row.append("\ttime_stamp")
            historyWriter.writerow(row)

    # --------------------------------------------------------------------------
    def LogCurrentValues( self, optimization_iteration, communicator, additional_values ):
        self.current_iteration = optimization_iteration
        if self.__IsFirstLog():
            self.initial_iteration = optimization_iteration
            self.__AddResponseValuesToHistory( communicator )
            self.__DetermineReferenceValuesForOutput( communicator )
        else:
            self.__AddResponseValuesToHistory( communicator )
            self.__DetermineChangeOfObjective()
        self.__PrintInfoAboutResponseFunctionValues()
        self.__WriteDataToLogFile( additional_values )
        self.previous_iteration = optimization_iteration

    # --------------------------------------------------------------------------
    def FinalizeLogging( self ):
        pass

    # --------------------------------------------------------------------------
    def GetValue( self ):
        pass

    # --------------------------------------------------------------------------
    def __CreateCompleteResponseLogFilename( self, optimization_settings ):
        resultsDirectory = optimization_settings["output"]["output_directory"].GetString()
        responseLogFilename = optimization_settings["output"]["response_log_filename"].GetString() + ".csv"
        return os.path.join( resultsDirectory, responseLogFilename )

    # -------------------------------------------------------------------------
    def __IsFirstLog( self ):
        if len(self.objective_history) == 0:
            return True
        else:
            return False

    # --------------------------------------------------------------------------
    def __AddResponseValuesToHistory( self, communicator ):
        J_id = self.specified_objectives[0]["identifier"].GetString()
        self.objective_history[self.current_iteration] = communicator.getValue( J_id )

        for itr in range(self.specified_constraints.size()):
            Ci_id = self.specified_constraints[itr]["identifier"].GetString()
            self.constraints_history[Ci_id][self.current_iteration] = communicator.getValue( Ci_id )

    # --------------------------------------------------------------------------
    def __DetermineReferenceValuesForOutput( self, communicator ):
        self.objective_output_reference = self.objective_history[self.initial_iteration]
        if abs(self.objective_output_reference)<1e-12:
            print("\n> WARNING: Objective reference value < 1e-12!! Therefore, standard reference value of 1 is assumed! ")
            self.objective_output_reference = 1.0
        else:
            self.objective_output_reference = self.objective_history[self.initial_iteration]

        for itr in range(self.specified_constraints.size()):
            Ci_id = self.specified_constraints[itr]["identifier"].GetString()
            self.constraints_output_reference[Ci_id] = communicator.getReferenceValue( Ci_id )

    # --------------------------------------------------------------------------
    def __DetermineChangeOfObjective( self ):
        objective_value = self.objective_history[self.current_iteration]
        previos_objective_value = self.objective_history[self.previous_iteration]

        self.abs_change_of_objective = 100*(objective_value-self.objective_output_reference) / abs(self.objective_output_reference)
        self.rel_change_of_objective = 100*(objective_value-previos_objective_value) / abs(self.objective_output_reference)

    # --------------------------------------------------------------------------
    def __PrintInfoAboutResponseFunctionValues( self ):
        print("\n> Current value of objective function = ", round(self.objective_history[self.current_iteration],12))
        print("> Absolut change of objective function = ",round(self.abs_change_of_objective,4)," [%]")
        print("> Relative change of objective function = ",round(self.rel_change_of_objective,4)," [%]\n")

        for itr in range(self.specified_constraints.size()):
            Ci_id = self.specified_constraints[itr]["identifier"].GetString()
            ci = self.constraints_history[Ci_id][self.current_iteration]
            print("> Value of C"+str(itr)+" = ", round(ci,12))

    # --------------------------------------------------------------------------
    def __WriteDataToLogFile( self, additional_values ):
        with open(self.complete_response_log_filename, 'a') as csvfile:
            historyWriter = csv.writer(csvfile, delimiter=',',quotechar='|',quoting=csv.QUOTE_MINIMAL)
            row = []
            row.append("\t"+str(self.current_iteration))
            row.append("\t"+str(self.objective_history[self.current_iteration]))
            row.append("\t"+str(self.abs_change_of_objective))
            row.append("\t"+str(self.rel_change_of_objective))

            for itr in range(self.specified_constraints.size()):
                Ci_id = self.specified_constraints[itr]["identifier"].GetString()
                ci = self.constraints_history[Ci_id][self.current_iteration]
                ci_ref = self.constraints_output_reference[Ci_id]
                row.append("\t"+str(ci))
                row.append("\t"+str(ci_ref))

            row.append("\t"+str(additional_values["step_length"]))
            row.append("\t"+Timer().GetTimeStamp())
            historyWriter.writerow(row)

# ==============================================================================
