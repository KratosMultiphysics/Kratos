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
from value_logger_base import ValueLogger

# Import additional libraries
import csv
from custom_timer import Timer

# ==============================================================================
class ValueLoggerTrustRegion( ValueLogger ):

    # --------------------------------------------------------------------------
    def __init__( self, optimization_settings ):
        self.optimization_settings = optimization_settings

        self.specified_objectives = optimization_settings["objectives"]
        self.specified_constraints = optimization_settings["constraints"]

        self.complete_log_file_name = self.__CreateCompleteLogFileName( optimization_settings )

        self.value_history = {}

        self.objective_reference_value = None

        self.current_iteration = 0
        self.previos_iteration = 0

    # --------------------------------------------------------------------------
    def InitializeLogging( self ):
        with open(self.complete_log_file_name, 'w') as csvfile:
            historyWriter = csv.writer(csvfile, delimiter=',',quotechar='|',quoting=csv.QUOTE_MINIMAL)
            row = []
            row.append("{:>4s}".format("itr"))
            row.append("{:>20s}".format("f"))
            row.append("{:>12s}".format("df_abs[%]"))
            row.append("{:>12s}".format("df_rel[%]"))

            for itr in range(self.specified_constraints.size()):
                con = self.specified_constraints[itr]
                con_type = con["type"].GetString()
                row.append("{:>20s}".format("c"+str(itr+1)+": "+con_type))
                row.append("{:>20s}".format("ref_c"+str(itr+1)))
                row.append("{:>12s}".format("len_c"+str(itr+1)))
                row.append("{:>12s}".format("adj_len_c"+str(itr+1)))

            row.append("{:>12s}".format("bi_itrs"))
            row.append("{:>12s}".format("bi_err"))
            row.append("{:>17}".format("test_norm_dx_bar"))
            row.append("{:>12s}".format("norm_dx"))
            row.append("{:>12s}".format("step_length"))
            row.append("{:>25s}".format("time_stamp"))
            historyWriter.writerow(row)

    # --------------------------------------------------------------------------
    def LogCurrentValues( self, optimization_iteration, communicator, additional_values ):
        self.current_iteration = optimization_iteration

        self.__AddValuesToHistory( communicator, additional_values )
        self.__PrintSelectedValuesToConsole( communicator, additional_values )
        self.__WriteDataToLogFile( communicator, additional_values )

        self.previos_iteration = optimization_iteration

    # --------------------------------------------------------------------------
    def FinalizeLogging( self ):
        pass

    # --------------------------------------------------------------------------
    def GetHistoryOfValues( self ):
        return self.value_history

    # --------------------------------------------------------------------------
    def __CreateCompleteLogFileName( self, optimization_settings ):
        resultsDirectory = optimization_settings["output"]["output_directory"].GetString()
        responseLogFilename = optimization_settings["output"]["response_log_filename"].GetString() + ".csv"
        return os.path.join( resultsDirectory, responseLogFilename )

    # --------------------------------------------------------------------------
    def __AddValuesToHistory( self, communicator, additional_values ):
        # Add objective values and change of objective values
        objective_id = self.specified_objectives[0]["identifier"].GetString()

        if self.__IsFirstLog(objective_id):
            self.value_history[objective_id] = {}
            self.value_history[objective_id][self.current_iteration] = communicator.getValue( objective_id )

            self.objective_reference_value = communicator.getValue( objective_id )
            if abs(self.objective_reference_value)<1e-12:
                print("\n> WARNING: Objective reference value < 1e-12!! Therefore, standard reference value of 1 is assumed! ")
                self.objective_reference_value = 1.0

            self.value_history["abs_change_obj"] = {self.current_iteration: 0.0}
            self.value_history["rel_change_obj"] = {self.current_iteration: 0.0}
        else:
            self.value_history[objective_id][self.current_iteration] = communicator.getValue( objective_id )

            current_obj_value = self.value_history[objective_id][self.current_iteration]
            previos_obj_value = self.value_history[objective_id][self.previos_iteration]
            abs_change = 100*(current_obj_value-self.objective_reference_value) / abs(self.objective_reference_value)
            rel_change = 100*(current_obj_value-previos_obj_value) / abs(self.objective_reference_value)

            self.value_history["abs_change_obj"][self.current_iteration] = abs_change
            self.value_history["rel_change_obj"][self.current_iteration] = rel_change

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

    # -------------------------------------------------------------------------
    def __IsFirstLog( self, key ):
        if key in self.value_history.keys():
            return False
        else:
            return True

    # --------------------------------------------------------------------------
    def __PrintSelectedValuesToConsole( self, communicator, additional_values ):
        print("\n-------------------------------------------------------")

        objective_id = self.specified_objectives[0]["identifier"].GetString()
        print("\n> Current value of objective = ", round(communicator.getValue(objective_id),12))

        print("> Absolut change of objective = ",round(self.value_history["abs_change_obj"][self.current_iteration],4)," [%]")
        print("> Relative change of objective = ",round(self.value_history["rel_change_obj"][self.current_iteration],4)," [%]\n")

        for itr in range(self.specified_constraints.size()):
            constraint_id = self.specified_constraints[itr]["identifier"].GetString()
            print("> Value of C"+str(itr)+" = ", round(communicator.getValue(constraint_id),12))

        print("\nNormInf3D of dx = ", round(additional_values["norm_dx"],6))

        print("\nlen_bar_obj = ", round(additional_values["len_bar_obj"],6))
        print("adj_len_bar_obj = ", round(additional_values["adj_len_bar_obj"],6))

        print("\nlen_bar_cons = ", [round(entry, 6) for entry in additional_values["len_bar_cons"]])
        print("adj_len_bar_cons = ", [round(entry, 6) for entry in additional_values["adj_len_bar_cons"]])

        print("\n-------------------------------------------------------")

    # --------------------------------------------------------------------------
    def __WriteDataToLogFile( self, communicator, additional_values ):

        with open(self.complete_log_file_name, 'a') as csvfile:
            historyWriter = csv.writer(csvfile, delimiter=',',quotechar='|',quoting=csv.QUOTE_MINIMAL)
            row = []
            row.append(str("{:>4d}".format(self.current_iteration)))

            objective_id = self.specified_objectives[0]["identifier"].GetString()
            row.append(str("{:>20f}".format(communicator.getValue(objective_id))))

            row.append(str("{:>12f}".format(self.value_history["abs_change_obj"][self.current_iteration])))
            row.append(str("{:>12f}".format(self.value_history["rel_change_obj"][self.current_iteration])))

            for itr in range(self.specified_constraints.size()):
                constraint_id = self.specified_constraints[itr]["identifier"].GetString()
                row.append(str("{:>20f}".format(communicator.getValue(constraint_id))))
                row.append(str("{:>20f}".format(communicator.getReferenceValue(constraint_id))))
                row.append(str("{:>12f}".format(additional_values["len_bar_cons"][itr])))
                row.append(str("{:>12f}".format(additional_values["adj_len_bar_cons"][itr])))

            row.append(str("{:>12d}".format(additional_values["bi_itrs"])))
            row.append(str("{:>12f}".format(additional_values["bi_err"])))
            row.append(str("{:>17f}".format(additional_values["test_norm_dx_bar"])))
            row.append(str("{:>12f}".format(additional_values["norm_dx"])))
            row.append(str("{:>12f}".format(additional_values["step_length"])))
            row.append("{:>25}".format(Timer().GetTimeStamp()))
            historyWriter.writerow(row)

# ==============================================================================
