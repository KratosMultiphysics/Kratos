# ==============================================================================
#  KratosShapeOptimizationApplication
#
#  License:         BSD License
#                   license: ShapeOptimizationApplication/license.txt
#
#  Main authors:    Baumgaertner Daniel, https://github.com/dbaumgaertner
#
# ==============================================================================

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
    def InitializeLogging( self ):
        with open(self.complete_log_file_name, 'w') as csvfile:
            historyWriter = csv.writer(csvfile, delimiter=',',quotechar='|',quoting=csv.QUOTE_MINIMAL)
            row = []
            row.append("{:>4s}".format("itr"))
            row.append("{:>13}".format("f"))
            row.append("{:>12s}".format("df_abs[%]"))
            row.append("{:>12s}".format("df_rel[%]"))

            for itr in range(self.specified_constraints.size()):
                con_type = self.specified_constraints[itr]["type"].GetString()
                row.append("{:>13}".format("c"+str(itr+1)+": "+con_type))
                row.append("{:>13}".format("c"+str(itr+1)+"_ref"))
                row.append("{:>12s}".format("len_C"+str(itr+1)))
                row.append("{:>12s}".format("adj_len_C"+str(itr+1)))

            row.append("{:>12s}".format("bi_itrs"))
            row.append("{:>12s}".format("bi_err"))
            row.append("{:>17}".format("test_norm_dX_bar"))
            row.append("{:>12s}".format("norm_dX"))
            row.append("{:>12s}".format("step_length"))
            row.append("{:>25s}".format("time_stamp"))
            historyWriter.writerow(row)

    # --------------------------------------------------------------------------
    def _WriteCurrentValuesToConsole( self ):
        print("\n-------------------------------------------------------")

        objective_id = self.specified_objectives[0]["identifier"].GetString()
        print("\n> Current value of objective = ", round(self.value_history[objective_id][self.current_iteration],12))

        print("> Absolut change of objective = ",round(self.value_history["abs_change_obj"][self.current_iteration],4)," [%]")
        print("> Relative change of objective = ",round(self.value_history["rel_change_obj"][self.current_iteration],4)," [%]\n")

        for itr in range(self.specified_constraints.size()):
            constraint_id = self.specified_constraints[itr]["identifier"].GetString()
            print("> Value of C"+str(itr+1)+" = ", round(self.value_history[constraint_id][self.current_iteration],12))

        print("\nNormInf3D of dX = ", round(self.value_history["norm_dX"][self.current_iteration],6))

        print("\nlen_bar_obj = ", round(self.value_history["len_bar_obj"][self.current_iteration],6))
        print("adj_len_bar_obj = ", round(self.value_history["adj_len_bar_obj"][self.current_iteration],6))

        print("\nlen_bar_cons = ", [round(entry, 6) for entry in self.value_history["len_bar_cons"][self.current_iteration]])
        print("adj_len_bar_cons = ", [round(entry, 6) for entry in self.value_history["adj_len_bar_cons"][self.current_iteration]])

        print("\n-------------------------------------------------------")

    # --------------------------------------------------------------------------
    def _WriteCurrentValuesToFile( self ):
        with open(self.complete_log_file_name, 'a') as csvfile:
            historyWriter = csv.writer(csvfile, delimiter=',',quotechar='|',quoting=csv.QUOTE_MINIMAL)
            row = []
            row.append("{:>4d}".format(self.current_iteration))

            objective_id = self.specified_objectives[0]["identifier"].GetString()
            row.append(" {:> .5E}".format(self.value_history[objective_id][self.current_iteration]))
            row.append("{:>12f}".format(self.value_history["abs_change_obj"][self.current_iteration]))
            row.append("{:>12f}".format(self.value_history["rel_change_obj"][self.current_iteration]))

            for itr in range(self.specified_constraints.size()):
                constraint_id = self.specified_constraints[itr]["identifier"].GetString()
                row.append(" {:> .5E}".format(self.value_history[constraint_id][self.current_iteration]))
                row.append(" {:> .5E}".format(self.communicator.getReferenceValue(constraint_id)))
                row.append("{:>12f}".format(self.value_history["len_bar_cons"][self.current_iteration][itr]))
                row.append("{:>12f}".format(self.value_history["adj_len_bar_cons"][self.current_iteration][itr]))

            row.append("{:>12d}".format(self.value_history["bi_itrs"][self.current_iteration]))
            row.append("{:>12f}".format(self.value_history["bi_err"][self.current_iteration]))
            row.append("{:>17f}".format(self.value_history["test_norm_dX_bar"][self.current_iteration]))
            row.append("{:>12f}".format(self.value_history["norm_dX"][self.current_iteration]))
            row.append(" {:>.5E}".format(self.value_history["step_length"][self.current_iteration]))
            row.append("{:>25}".format(Timer().GetTimeStamp()))
            historyWriter.writerow(row)

# ==============================================================================
