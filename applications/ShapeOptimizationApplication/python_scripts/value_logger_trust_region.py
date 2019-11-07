# ==============================================================================
#  KratosShapeOptimizationApplication
#
#  License:         BSD License
#                   license: ShapeOptimizationApplication/license.txt
#
#  Main authors:    Baumgaertner Daniel, https://github.com/dbaumgaertner
#
# ==============================================================================

# importing the Kratos Library
import KratosMultiphysics as KM

# Import logger base classes
from .value_logger_base import ValueLogger

# Import additional libraries
import csv
from .custom_timer import Timer

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

            for itr in range(self.constraints.size()):
                con_type = self.constraints[itr]["type"].GetString()
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
        KM.Logger.Print("")
        KM.Logger.Print("-------------------------------------------------------\n")

        objective_id = self.objectives[0]["identifier"].GetString()
        KM.Logger.PrintInfo("ShapeOpt", "Current value of objective = ", round(self.history["response_value"][objective_id][self.current_index],12))

        KM.Logger.PrintInfo("ShapeOpt", "Absolut change of objective = ",round(self.history["abs_change_objective"][self.current_index],4)," [%]")
        KM.Logger.PrintInfo("ShapeOpt", "Relative change of objective = ",round(self.history["rel_change_objective"][self.current_index],4)," [%]\n")

        for itr in range(self.constraints.size()):
            constraint_id = self.constraints[itr]["identifier"].GetString()
            KM.Logger.PrintInfo("ShapeOpt", "Value of C"+str(itr+1)+" = ", round(self.history["response_value"][constraint_id][self.current_index],12), "\n")

        KM.Logger.PrintInfo("ShapeOpt", "NormInf3D of dX = ", round(self.history["norm_dX"][self.current_index],6), "\n")

        KM.Logger.PrintInfo("ShapeOpt", "len_bar_obj = ", round(self.history["len_bar_obj"][self.current_index],6))
        KM.Logger.PrintInfo("ShapeOpt", "adj_len_bar_obj = ", round(self.history["adj_len_bar_obj"][self.current_index],6), "\n")

        KM.Logger.PrintInfo("ShapeOpt", "len_bar_cons = ", [round(entry, 6) for entry in self.history["len_bar_cons"][self.current_index]])
        KM.Logger.PrintInfo("ShapeOpt", "adj_len_bar_cons = ", [round(entry, 6) for entry in self.history["adj_len_bar_cons"][self.current_index]], "\n")

        KM.Logger.Print("-------------------------------------------------------")

    # --------------------------------------------------------------------------
    def _WriteCurrentValuesToFile( self ):
        with open(self.complete_log_file_name, 'a') as csvfile:
            historyWriter = csv.writer(csvfile, delimiter=',',quotechar='|',quoting=csv.QUOTE_MINIMAL)
            row = []
            row.append("{:>4d}".format(self.current_index))

            objective_id = self.objectives[0]["identifier"].GetString()
            row.append(" {:> .5E}".format(self.history["response_value"][objective_id][self.current_index]))
            row.append(" {:> .5E}".format(self.history["abs_change_objective"][self.current_index]))
            row.append(" {:> .5E}".format(self.history["rel_change_objective"][self.current_index]))

            for itr in range(self.constraints.size()):
                constraint_id = self.constraints[itr]["identifier"].GetString()
                row.append(" {:> .5E}".format(self.history["response_value"][constraint_id][self.current_index]))
                row.append(" {:> .5E}".format(self.communicator.getReferenceValue(constraint_id)))
                row.append("{:>12f}".format(self.history["len_bar_cons"][self.current_index][itr]))
                row.append("{:>12f}".format(self.history["adj_len_bar_cons"][self.current_index][itr]))

            row.append("{:>12d}".format(self.history["bi_itrs"][self.current_index]))
            row.append("{:>12f}".format(self.history["bi_err"][self.current_index]))
            row.append("{:>17f}".format(self.history["test_norm_dX_bar"][self.current_index]))
            row.append("{:>12f}".format(self.history["norm_dX"][self.current_index]))
            row.append(" {:>.5E}".format(self.history["step_length"][self.current_index]))
            row.append("{:>25}".format(Timer().GetTimeStamp()))
            historyWriter.writerow(row)

# ==============================================================================
