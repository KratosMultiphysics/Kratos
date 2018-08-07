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

        self.objective_id = optimization_settings["objectives"][0]["identifier"].GetString()

        self.complete_response_log_filename = self.__CreateCompleteResponseLogFilename( optimization_settings )

        self.value_history = {}

        self.objective_reference_value = None
        self.constraints_reference_values = {}

        self.abs_change_of_objective = 0.0
        self.rel_change_of_objective = 0.0

        self.current_iteration = 0
        self.previos_iteration = 0

    # --------------------------------------------------------------------------
    def InitializeLogging( self ):
        self.__InitializeResponseFile()

    # --------------------------------------------------------------------------
    def LogCurrentValues( self, optimization_iteration, communicator, additional_values ):
        self.current_iteration = optimization_iteration

        if self.__IsFirstLog( self.objective_id ):
            self.__AddValuesToHistory( communicator, additional_values )
            self.__DetermineReferenceValuesForOutput( communicator )
        else:
            self.__AddValuesToHistory( communicator, additional_values )
            self.__DetermineChangeOfObjective()
        self.__PrintSelectedValuesToConsole()
        self.__WriteDataToLogFile( communicator, additional_values )

        self.previos_iteration = optimization_iteration

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

    # --------------------------------------------------------------------------
    def __InitializeResponseFile( self ):
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

            row.append("\tnorm_dx[-]")
            row.append("\tstep_length[-]")
            row.append("\ttime_stamp")
            historyWriter.writerow(row)

    # -------------------------------------------------------------------------
    def __IsFirstLog( self, key ):
        if key in self.value_history.keys():
            return False
        else:
            return True

    # --------------------------------------------------------------------------
    def __AddValuesToHistory( self, communicator, additional_values ):
        # Add objective values
        if self.__IsFirstLog(self.objective_id):
            self.value_history[self.objective_id] = {}
        self.value_history[self.objective_id][self.current_iteration] = communicator.getValue( self.objective_id )

        # Add constraint values
        for itr in range(self.specified_constraints.size()):
            constraint_id = self.specified_constraints[itr]["identifier"].GetString()
            if self.__IsFirstLog(constraint_id):
                self.value_history[constraint_id] = {}
            self.value_history[constraint_id][self.current_iteration] = communicator.getValue( constraint_id )

        # Add additional values
        for key, value in additional_values.items():
            if self.__IsFirstLog(key):
                self.value_history[key] = {}
            self.value_history[key][self.current_iteration] = value

    # --------------------------------------------------------------------------
    def __DetermineReferenceValuesForOutput( self, communicator ):
        self.objective_reference_value = communicator.getValue( self.objective_id )

        if abs(self.objective_reference_value)<1e-12:
            print("\n> WARNING: Objective reference value < 1e-12!! Therefore, standard reference value of 1 is assumed! ")
            self.objective_reference_value = 1.0

        for itr in range(self.specified_constraints.size()):
            constraint_id = self.specified_constraints[itr]["identifier"].GetString()
            self.constraints_reference_values[constraint_id] = communicator.getReferenceValue( constraint_id )

    # --------------------------------------------------------------------------
    def __DetermineChangeOfObjective( self ):
        objective_value = self.value_history[self.objective_id][self.current_iteration]
        previos_objective_value = self.value_history[self.objective_id][self.previos_iteration]

        self.abs_change_of_objective = 100*(objective_value-self.objective_reference_value) / abs(self.objective_reference_value)
        self.rel_change_of_objective = 100*(objective_value-previos_objective_value) / abs(self.objective_reference_value)

    # --------------------------------------------------------------------------
    def __PrintSelectedValuesToConsole( self ):
        print("\n-------------------------------------------------------")

        print("\n> Current value of objective = ", round(self.value_history[self.objective_id][self.current_iteration],12))
        print("> Absolut change of objective = ",round(self.abs_change_of_objective,4)," [%]")
        print("> Relative change of objective = ",round(self.rel_change_of_objective,4)," [%]\n")

        for itr in range(self.specified_constraints.size()):
            Ci_id = self.specified_constraints[itr]["identifier"].GetString()
            ci = self.value_history[Ci_id][self.current_iteration]
            print("> Value of C"+str(itr)+" = ", round(ci,12))

        print("\nNormInf3D of dx = ", self.value_history["norm_dx"][self.current_iteration])

        print("\nlen_bar_obj = ", self.value_history["len_bar_obj"][self.current_iteration])
        print("adj_len_bar_obj = ", self.value_history["adj_len_bar_obj"][self.current_iteration])

        print("\nlen_bar_ineqs = ", self.value_history["len_bar_ineqs"][self.current_iteration])
        print("adj_len_bar_ineqs = ", self.value_history["adj_len_bar_ineqs"][self.current_iteration])

        print("\nlen_bar_eqs = ", self.value_history["len_bar_eqs"][self.current_iteration])
        print("adj_len_bar_eqs = ", self.value_history["adj_len_bar_eqs"][self.current_iteration])

        print("\n-------------------------------------------------------")

    # --------------------------------------------------------------------------
    def __WriteDataToLogFile( self, communicator, additional_values ):
        with open(self.complete_response_log_filename, 'a') as csvfile:
            historyWriter = csv.writer(csvfile, delimiter=',',quotechar='|',quoting=csv.QUOTE_MINIMAL)
            row = []
            row.append("\t"+str(self.current_iteration))
            row.append("\t"+str(self.value_history[self.objective_id][self.current_iteration]))
            row.append("\t"+str(self.abs_change_of_objective))
            row.append("\t"+str(self.rel_change_of_objective))

            for itr in range(self.specified_constraints.size()):
                Ci_id = self.specified_constraints[itr]["identifier"].GetString()
                ci = self.value_history[Ci_id][self.current_iteration]
                ci_ref = self.constraints_reference_values[Ci_id]
                row.append("\t"+str(ci))
                row.append("\t"+str(ci_ref))

            row.append("\t"+str(additional_values["norm_dx"]))
            row.append("\t"+str(additional_values["step_length"]))
            row.append("\t"+Timer().GetTimeStamp())
            historyWriter.writerow(row)

# ==============================================================================
