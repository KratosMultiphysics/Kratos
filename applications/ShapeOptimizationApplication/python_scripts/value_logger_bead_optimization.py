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

# Kratos Core and Apps
from KratosMultiphysics import *
from KratosMultiphysics.ShapeOptimizationApplication import *

# Import logger base classes
from value_logger_base import ValueLogger

# Import additional libraries
import csv
from custom_timer import Timer

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
            row.append("{:>12s}".format("dl_rel[%]"))
            row.append("{:>13s}".format("f"))
            row.append("{:>12s}".format("df_abs[%]"))
            row.append("{:>12s}".format("df_rel[%]"))
            row.append("{:>13s}".format("lambda_p"))
            row.append("{:>13s}".format("p"))
            row.append("{:>13s}".format("p_scaling"))
            row.append("{:>13s}".format("p_fac"))
            row.append("{:>12s}".format("step_size"))
            row.append("{:>25s}".format("time_stamp"))
            historyWriter.writerow(row)

    # --------------------------------------------------------------------------
    def _WriteCurrentValuesToConsole( self ):
        objective_id = self.specified_objectives[0]["identifier"].GetString()
        print("\n> Current value of objective = ", round(self.value_history[objective_id][self.current_iteration],12))

        print("> Absolut change of objective = ",round(self.value_history["abs_change_obj"][self.current_iteration],4)," [%]")
        print("> Relative change of objective = ",round(self.value_history["rel_change_obj"][self.current_iteration],4)," [%]\n")
        print("> Value of penalty term = ",round(self.value_history["penalty_value"][self.current_iteration],12),"\n")

    # --------------------------------------------------------------------------
    def _WriteCurrentValuesToFile( self ):
        with open(self.complete_log_file_name, 'a') as csvfile:
            historyWriter = csv.writer(csvfile, delimiter=',',quotechar='|',quoting=csv.QUOTE_MINIMAL)
            row = []
            row.append("{:>12d}".format(self.current_iteration))
            row.append("{:>10d}".format(self.value_history["outer_iteration"][self.current_iteration]))
            row.append("{:>10d}".format(self.value_history["inner_iteration"][self.current_iteration]))
            objective_id = self.specified_objectives[0]["identifier"].GetString()
            row.append(" {:> .5E}".format(self.value_history["lagrange_value"][self.current_iteration]))
            row.append("{:>12f}".format(self.value_history["lagrange_value_relative_change"][self.current_iteration]))
            row.append(" {:> .5E}".format(self.value_history[objective_id][self.current_iteration]))
            row.append("{:>12f}".format(self.value_history["abs_change_obj"][self.current_iteration]))
            row.append("{:>12f}".format(self.value_history["rel_change_obj"][self.current_iteration]))
            row.append(" {:> .5E}".format(self.value_history["penalty_lambda"][self.current_iteration]))
            row.append(" {:> .5E}".format(self.value_history["penalty_value"][self.current_iteration]))
            row.append(" {:> .5E}".format(self.value_history["penalty_scaling"][self.current_iteration]))
            row.append(" {:> .5E}".format(self.value_history["penalty_factor"][self.current_iteration]))
            row.append("{:>12f}".format(self.value_history["step_size"][self.current_iteration]))
            row.append("{:>25}".format(Timer().GetTimeStamp()))
            historyWriter.writerow(row)

# ==============================================================================
