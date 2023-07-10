# ==============================================================================
#  KratosShapeOptimizationApplication
#
#  License:         BSD License
#                   license: ShapeOptimizationApplication/license.txt
#
#  Main authors:    Ihar Antonau
#
# ==============================================================================

# importing the Kratos Library
import KratosMultiphysics as Kratos

# Import logger base classes
from .value_logger_base import ValueLogger

# Import additional libraries
import csv
from KratosMultiphysics.ShapeOptimizationApplication.utilities.custom_timer import Timer

# ==============================================================================
class ValueLoggerRelaxedGradientProjection( ValueLogger ):
    # --------------------------------------------------------------------------
    def InitializeLogging( self ):
        with open(self.complete_log_file_name, 'w') as csvfile:
            historyWriter = csv.writer(csvfile, delimiter=',',quotechar='|',quoting=csv.QUOTE_MINIMAL)
            row = []
            row.append("{:>4s}".format("itr"))
            row.append("{:>13s}".format("f"))
            row.append("{:>13s}".format("df_abs[%]"))
            row.append("{:>13s}".format("df_rel[%]"))

            for itr in range(self.constraints.size()):
                con_type = self.constraints[itr]["type"].GetString()
                row.append("{:>13s}".format("c"+str(itr+1)+": "+con_type))
                row.append("{:>13s}".format("c"+str(itr+1)+"_ref"))

            row.append("{:>13s}".format("step_size"))
            row.append("{:>13s}".format("inf_norm_p"))
            row.append("{:>13s}".format("inf_norm_c"))
            row.append("{:>13s}".format("projection_norm"))

            for itr in range(self.constraints.size()):
                row.append("{:>13s}".format("c"+str(itr+1)+"_buffer_value"))
                row.append("{:>13s}".format("c"+str(itr+1)+"_buffer_size"))
                row.append("{:>13s}".format("c"+str(itr+1)+"_buffer_size_factor"))
                row.append("{:>13s}".format("c"+str(itr+1)+"_central_buffer_value"))
                row.append("{:>13s}".format("c"+str(itr+1)+"_lower_buffer_value"))
                row.append("{:>13s}".format("c"+str(itr+1)+"_upper_buffer_value"))

            row.append("{:>25s}".format("time_stamp"))
            historyWriter.writerow(row)

    # --------------------------------------------------------------------------
    def _WriteCurrentValuesToConsole( self ):
        objective_id = self.objectives[0]["identifier"].GetString()
        Kratos.Logger.Print("")
        Kratos.Logger.PrintInfo("ShapeOpt", "Current value of objective = ", "{:> .5E}".format(self.history["response_value"]
                                                                                               [objective_id][self.current_index]))

        Kratos.Logger.PrintInfo("ShapeOpt", "Absolute change of objective = ","{:> .5E}".format(self.history["abs_change_objective"][self.current_index]),
                                                                                                " [%]")
        Kratos.Logger.PrintInfo("ShapeOpt", "Relative change of objective = ","{:> .5E}".format(self.history["rel_change_objective"][self.current_index]),
                                                                                                " [%]\n")

        for itr in range(self.constraints.size()):
            constraint_id = self.constraints[itr]["identifier"].GetString()
            Kratos.Logger.PrintInfo("ShapeOpt", f"Value of C{str(itr+1)} = ","{:> .5E}".format(
                                    self.history["response_value"][constraint_id][self.current_index]))

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

            row.append(" {:> .5E}".format(self.history["step_size"][self.current_index]))
            row.append(" {:> .5E}".format(self.history["inf_norm_p"][self.current_index]))
            row.append(" {:> .5E}".format(self.history["inf_norm_c"][self.current_index]))
            row.append(" {:> .5E}".format(self.history["projection_norm"][self.current_index]))

            for itr in range(self.constraints.size()):
                constraint_id = self.constraints[itr]["identifier"].GetString()
                row.append(" {:> .5E}".format(self.history["c"+str(itr+1)+"_buffer_value"][self.current_index]))
                row.append(" {:> .5E}".format(self.history["c"+str(itr+1)+"_buffer_size"][self.current_index]))
                row.append(" {:> .5E}".format(self.history["c"+str(itr+1)+"_buffer_size_factor"][self.current_index]))
                row.append(" {:> .5E}".format(self.history["c"+str(itr+1)+"_central_buffer_value"][self.current_index] +
                                              self.communicator.getReferenceValue(constraint_id)))
                row.append(" {:> .5E}".format(self.history["c"+str(itr+1)+"_lower_buffer_value"][self.current_index]
                                              + self.communicator.getReferenceValue(constraint_id)))
                row.append(" {:> .5E}".format(self.history["c"+str(itr+1)+"_upper_buffer_value"][self.current_index]
                                              + self.communicator.getReferenceValue(constraint_id)))

            row.append("{:>25}".format(Timer().GetTimeStamp()))
            historyWriter.writerow(row)

# ==============================================================================