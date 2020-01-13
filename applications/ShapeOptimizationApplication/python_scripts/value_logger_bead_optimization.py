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

# importing the Kratos Library
import KratosMultiphysics as KM

# Import logger base classes
from .value_logger_base import ValueLogger

# Import additional libraries
import csv
from .custom_timer import Timer

# ==============================================================================
class ValueLoggerBeadOptimization( ValueLogger ):
    # --------------------------------------------------------------------------
    def InitializeLogging( self ):
        with open(self.complete_log_file_name, 'w') as csvfile:
            historyWriter = csv.writer(csvfile, delimiter=',',quotechar='|',quoting=csv.QUOTE_MINIMAL)
            row = []
            row.append("{:>12s}".format("overall_itr"))
            row.append("{:>10s}".format("outer_itr"))
            row.append("{:>10s}".format("inner_itr"))
            row.append("{:>13s}".format("l"))
            row.append("{:>13s}".format("dl_rel[%]"))
            row.append("{:>13s}".format("f"))
            row.append("{:>13s}".format("df_abs[%]"))
            row.append("{:>13s}".format("df_rel[%]"))
            row.append("{:>14s}".format("max_norm_dfda"))
            row.append("{:>13s}".format("lambda_p"))
            row.append("{:>13s}".format("p"))
            row.append("{:>13s}".format("p_scaling"))
            row.append("{:>13s}".format("p_fac"))
            row.append("{:>13s}".format("step_size"))
            row.append("{:>25s}".format("time_stamp"))
            historyWriter.writerow(row)

    # --------------------------------------------------------------------------
    def _WriteCurrentValuesToConsole( self ):
        objective_id = self.objectives[0]["identifier"].GetString()
        KM.Logger.Print("")
        KM.Logger.PrintInfo("ShapeOpt", "Current value of objective = ","{:> .5E}".format(self.history["response_value"][objective_id][self.current_index]))

        KM.Logger.PrintInfo("ShapeOpt", "Absolut change of objective = ","{:> .5E}".format(self.history["abs_change_objective"][self.current_index])," [%]")
        KM.Logger.PrintInfo("ShapeOpt", "Relative change of objective = ","{:> .5E}".format(self.history["rel_change_objective"][self.current_index])," [%]\n")
        KM.Logger.PrintInfo("ShapeOpt", "Value of penalty term = ","{:> .5E}".format(self.history["penalty_value"][self.current_index]),"\n")

    # --------------------------------------------------------------------------
    def _WriteCurrentValuesToFile( self ):
        with open(self.complete_log_file_name, 'a') as csvfile:
            historyWriter = csv.writer(csvfile, delimiter=',',quotechar='|',quoting=csv.QUOTE_MINIMAL)
            row = []
            row.append("{:>12d}".format(self.current_index))
            row.append("{:>10d}".format(self.history["outer_iteration"][self.current_index]))
            row.append("{:>10d}".format(self.history["inner_iteration"][self.current_index]))
            objective_id = self.objectives[0]["identifier"].GetString()
            row.append(" {:> .5E}".format(self.history["lagrange_value"][self.current_index]))
            row.append(" {:> .5E}".format(self.history["lagrange_value_relative_change"][self.current_index]))
            row.append(" {:> .5E}".format(self.history["response_value"][objective_id][self.current_index]))
            row.append(" {:> .5E}".format(self.history["abs_change_objective"][self.current_index]))
            row.append(" {:> .5E}".format(self.history["rel_change_objective"][self.current_index]))
            row.append("  {:> .5E}".format(self.history["max_norm_objective_gradient"][self.current_index]))
            row.append(" {:> .5E}".format(self.history["penalty_lambda"][self.current_index]))
            row.append(" {:> .5E}".format(self.history["penalty_value"][self.current_index]))
            row.append(" {:> .5E}".format(self.history["penalty_scaling"][self.current_index]))
            row.append(" {:> .5E}".format(self.history["penalty_factor"][self.current_index]))
            row.append(" {:> .5E}".format(self.history["step_size"][self.current_index]))
            row.append("{:>25}".format(Timer().GetTimeStamp()))
            historyWriter.writerow(row)

# ==============================================================================
