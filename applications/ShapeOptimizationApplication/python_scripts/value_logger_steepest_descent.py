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
class ValueLoggerSteepestDescent( ValueLogger ):
    # --------------------------------------------------------------------------
    def InitializeLogging( self ):
        with open(self.complete_log_file_name, 'w') as csvfile:
            historyWriter = csv.writer(csvfile, delimiter=',',quotechar='|',quoting=csv.QUOTE_MINIMAL)
            row = []
            row.append("{:>4s}".format("itr"))
            row.append("{:>13s}".format("f"))
            row.append("{:>13s}".format("df_abs[%]"))
            row.append("{:>13s}".format("df_rel[%]"))
            row.append("{:>13s}".format("step_size"))
            row.append("{:>25s}".format("time_stamp"))
            historyWriter.writerow(row)

    # --------------------------------------------------------------------------
    def _WriteCurrentValuesToConsole( self ):
        objective_id = self.specified_objectives[0]["identifier"].GetString()
        print("\n> Current value of objective = ", "{:> .5E}".format(self.value_history[objective_id][self.current_iteration]))

        print("> Absolut change of objective = ","{:> .5E}".format(self.value_history["abs_change_obj"][self.current_iteration])," [%]")
        print("> Relative change of objective = ","{:> .5E}".format(self.value_history["rel_change_obj"][self.current_iteration])," [%]\n")

    # --------------------------------------------------------------------------
    def _WriteCurrentValuesToFile( self ):
        with open(self.complete_log_file_name, 'a') as csvfile:
            historyWriter = csv.writer(csvfile, delimiter=',',quotechar='|',quoting=csv.QUOTE_MINIMAL)
            row = []
            row.append("{:>4d}".format(self.current_iteration))

            objective_id = self.specified_objectives[0]["identifier"].GetString()
            row.append(" {:> .5E}".format(self.value_history[objective_id][self.current_iteration]))
            row.append(" {:> .5E}".format(self.value_history["abs_change_obj"][self.current_iteration]))
            row.append(" {:> .5E}".format(self.value_history["rel_change_obj"][self.current_iteration]))

            row.append(" {:> .5E}".format(self.value_history["step_size"][self.current_iteration]))
            row.append("{:>25}".format(Timer().GetTimeStamp()))
            historyWriter.writerow(row)

# ==============================================================================
