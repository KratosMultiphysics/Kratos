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

# Import additional libraries
import os

# ==============================================================================
class ValueLogger():

    # --------------------------------------------------------------------------
    def __init__( self, communicator, optimization_settings ):
        self.communicator = communicator
        self.optimization_settings = optimization_settings

        self.objectives = optimization_settings["objectives"]
        self.constraints = optimization_settings["constraints"]

        self.complete_log_file_name = self.__CreateCompleteLogFileName( optimization_settings )

        self.obj_reference_value = 0
        self.current_index = 0
        self.previos_index = 0

        self.history = { "response_value"              : {},
                         "standardized_response_value" : {},
                         "abs_change_objective"        : {},
                         "rel_change_objective"        : {} }

        for itr in range(self.objectives.size()):
            objective_id = self.objectives[itr]["identifier"].GetString()
            self.history["response_value"][objective_id] = {}
            self.history["standardized_response_value"][objective_id] = {}

        for itr in range(self.constraints.size()):
            constraint_id = self.constraints[itr]["identifier"].GetString()
            self.history["response_value"][constraint_id] = {}
            self.history["standardized_response_value"][constraint_id] = {}

        self.predefined_keys = list(self.history.keys())

    # --------------------------------------------------------------------------
    def InitializeLogging( self ):
        pass

    # --------------------------------------------------------------------------
    def LogCurrentValues( self, index, additional_values_dictionary ):
        self.current_index = index

        self.__LogObjectiveValuesToHistory()
        self.__LogConstraintValuesToHistory()
        self.__LogAdditionalValuesToHistory( additional_values_dictionary )

        self._WriteCurrentValuesToConsole()
        self._WriteCurrentValuesToFile()

        self.previos_index = index

    # --------------------------------------------------------------------------
    def FinalizeLogging( self ):
        pass

    # --------------------------------------------------------------------------
    def GetValues( self, key ):
        return self.history[key]

    # --------------------------------------------------------------------------
    @staticmethod
    def __CreateCompleteLogFileName( optimization_settings ):
        resultsDirectory = optimization_settings["output"]["output_directory"].GetString()
        optimizationLogFilename = optimization_settings["output"]["optimization_log_filename"].GetString() + ".csv"
        return os.path.join( resultsDirectory, optimizationLogFilename )

    # --------------------------------------------------------------------------
    def __LogObjectiveValuesToHistory( self ):
        objective_id = self.objectives[0]["identifier"].GetString()
        is_first_log = self.history["response_value"][objective_id] == {}

        self.history["response_value"][objective_id][self.current_index] = self.communicator.getValue( objective_id )
        self.history["standardized_response_value"][objective_id][self.current_index] = self.communicator.getStandardizedValue( objective_id )

        if is_first_log:
            self.obj_reference_value = self.communicator.getValue( objective_id )

            if abs(self.obj_reference_value)<1e-12:
                KM.Logger.PrintWarning("ShapeOpt::ValueLogger", "Objective reference value < 1e-12!! Therefore, standard reference value of 1 is assumed! ")
                self.obj_reference_value = 1.0

            self.history["abs_change_objective"] = {self.current_index: 0.0}
            self.history["rel_change_objective"] = {self.current_index: 0.0}
        else:
            current_obj_value = self.history["response_value"][objective_id][self.current_index]
            previos_obj_value = self.history["response_value"][objective_id][self.previos_index]
            abs_change = 100*(current_obj_value-self.obj_reference_value) / abs(self.obj_reference_value)
            rel_change = 100*(current_obj_value-previos_obj_value) / abs(self.obj_reference_value)

            self.history["abs_change_objective"][self.current_index] = abs_change
            self.history["rel_change_objective"][self.current_index] = rel_change

    # --------------------------------------------------------------------------
    def __LogConstraintValuesToHistory( self ):
        for itr in range(self.constraints.size()):
            constraint_id = self.constraints[itr]["identifier"].GetString()
            self.history["response_value"][constraint_id][self.current_index] = self.communicator.getValue( constraint_id )
            self.history["standardized_response_value"][constraint_id][self.current_index] = self.communicator.getStandardizedValue( constraint_id )

    # --------------------------------------------------------------------------
    def __LogAdditionalValuesToHistory( self, additional_values_dictionary ):
        for key, value in additional_values_dictionary.items():
            if key in self.predefined_keys:
                raise NameError("ValueLogger: The key \""+key+"\" is already predefined and may not be overwritten! Predefined keys are: "+str(self.predefined_keys))

            if key not in self.history.keys():
                self.history[key] = {}

            self.history[key][self.current_index] = value

# ==============================================================================
