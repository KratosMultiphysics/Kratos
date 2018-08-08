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

# Import additional libraries
import os

# ==============================================================================
class ValueLogger():

    # --------------------------------------------------------------------------
    def __init__( self, communicator, optimization_settings ):
        self.communicator = communicator
        self.optimization_settings = optimization_settings

        self.specified_objectives = optimization_settings["objectives"]
        self.specified_constraints = optimization_settings["constraints"]

        self.complete_log_file_name = self.__CreateCompleteLogFileName( optimization_settings )

        self.obj_reference_value = 0
        self.current_iteration = 0
        self.previos_iteration = 0

        self.value_history = {}

    # --------------------------------------------------------------------------
    def InitializeLogging( self ):
        pass

    # --------------------------------------------------------------------------
    def AddValuesToHistory( self, current_iteration, additional_values ):
        self.current_iteration = current_iteration

        self.__AddObjectiveValuesToHistory()
        self.__AddConstraintValuesToHistory()
        self.__AddAdditionalValuesToHistory( additional_values )

        self.previos_iteration = current_iteration

    # --------------------------------------------------------------------------
    def LogCurrentValuesToConsole( self ):
        pass

    # --------------------------------------------------------------------------
    def LogCurrentValuesToFile( self ):
        pass

    # --------------------------------------------------------------------------
    def FinalizeLogging( self ):
        pass

    # --------------------------------------------------------------------------
    def GetHistoryOfLoggedValues( self ):
        return self.value_history

    # --------------------------------------------------------------------------
    def __CreateCompleteLogFileName( self, optimization_settings ):
        resultsDirectory = optimization_settings["output"]["output_directory"].GetString()
        responseLogFilename = optimization_settings["output"]["response_log_filename"].GetString() + ".csv"
        return os.path.join( resultsDirectory, responseLogFilename )

    # --------------------------------------------------------------------------
    def __AddObjectiveValuesToHistory( self ):
        objective_id = self.specified_objectives[0]["identifier"].GetString()

        if self.__IsFirstLog(objective_id):
            self.value_history[objective_id] = {}
            self.value_history[objective_id][self.current_iteration] = self.communicator.getValue( objective_id )

            self.obj_reference_value = self.communicator.getValue( objective_id )
            if abs(self.obj_reference_value)<1e-12:
                print("\n> WARNING: Objective reference value < 1e-12!! Therefore, standard reference value of 1 is assumed! ")
                self.obj_reference_value = 1.0

            self.value_history["abs_change_obj"] = {self.current_iteration: 0.0}
            self.value_history["rel_change_obj"] = {self.current_iteration: 0.0}
        else:
            self.value_history[objective_id][self.current_iteration] = self.communicator.getValue( objective_id )

            current_obj_value = self.value_history[objective_id][self.current_iteration]
            previos_obj_value = self.value_history[objective_id][self.previos_iteration]
            abs_change = 100*(current_obj_value-self.obj_reference_value) / abs(self.obj_reference_value)
            rel_change = 100*(current_obj_value-previos_obj_value) / abs(self.obj_reference_value)

            self.value_history["abs_change_obj"][self.current_iteration] = abs_change
            self.value_history["rel_change_obj"][self.current_iteration] = rel_change

    # --------------------------------------------------------------------------
    def __AddConstraintValuesToHistory( self ):
        for itr in range(self.specified_constraints.size()):
            constraint_id = self.specified_constraints[itr]["identifier"].GetString()

            if self.__IsFirstLog(constraint_id):
                self.value_history[constraint_id] = {}
            self.value_history[constraint_id][self.current_iteration] = self.communicator.getValue( constraint_id )

    # --------------------------------------------------------------------------
    def __AddAdditionalValuesToHistory( self, additional_values ):
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

# ==============================================================================
