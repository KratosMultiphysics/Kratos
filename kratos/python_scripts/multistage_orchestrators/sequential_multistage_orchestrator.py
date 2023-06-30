import KratosMultiphysics
from KratosMultiphysics.multistage_orchestrators.multistage_orchestrator import MultistageOrchestrator

class SequentialMultistageOrchestrator(MultistageOrchestrator):

    def __init__(self, settings) -> None:
        super().__init__(settings)

    def Run(self):
        '''Main function that runs a sequential multistage simulation.'''

        # Set the stages to be saved first
        # We do it first to warn the user about the non-present stages before running the first one
        self.__GetStagesToCheckpointList()

        # Check if loading from checkpoint is required and load
        if self.settings["orchestrator"]["settings"].Has("load_from_checkpoint"):
            if not self.settings["orchestrator"]["settings"]["load_from_checkpoint"].IsNull():
                if not self.settings["orchestrator"]["settings"]["load_from_checkpoint"].IsString():
                    err_msg = "'load_from_checkpoint' is expected to be null or a string with the loading point."
                    raise TypeError(err_msg)
                else:
                    loading_point = self.settings["orchestrator"]["settings"]["load_from_checkpoint"].GetString()
                    self.GetProject().Load(loading_point)

        # Run the stages list
        for stage_name in self.__GetExecutionList():
            # Check current stage settings
            self.CheckStageSettings(stage_name)

            # Set up and check current stage instance
            current_stage = self.CreateStage(stage_name)

            # Execute current stage preprocess
            self.RunCurrentStagePreprocess(stage_name)

            # Run (solve) current stage
            current_stage.Run()

            # Get the final data dictionary
            self.GetProject().output_data[stage_name] = current_stage.GetAnalysisStageFinalData()

            # Execute current stage postprocess
            self.RunCurrentStagePostprocess(stage_name)

            # Delete current stage instance
            del current_stage

            # Check the current stage is to be checkpointed and save it if so
            if stage_name in self.__GetStagesToCheckpointList():
                save_folder_name = self.__GetStagesCheckpointFolder()
                self.GetProject().Save(save_folder_name, stage_name)
                
            # Output validated simulation settings (with defaults) as simulation report
            self.__OutputValidatedSettings()

    def GetNumberOfStages(self):
        '''Returns the number of stages.'''

        return len(self.__GetExecutionList())

    def __GetExecutionList(self):
        '''Creates and returns the execution list.
        This method creates the execution list, either from a user-defined execution list or,
        if this is not provided, from the stages declaration order.
        '''

        # Check if the execution list has been already created
        if not hasattr(self, "__execution_list"):
            # Default case in which the execution list is provided by the user
            if self.settings["orchestrator"]["settings"].Has("execution_list"):
                self.__execution_list = self.settings["orchestrator"]["settings"]["execution_list"].GetStringArray()
            # If not provided, create an auxiliary execution list from the stages declaration order
            else:
                KratosMultiphysics.Logger.PrintInfo("'execution_list' is not provided. Stages will be executed according to their declaration order.")
                self.__execution_list = list(self.settings["stages"].keys())

        return self.__execution_list

    def __GetStagesToCheckpointList(self):
        '''Creates and returns the stages checkpoint list.

        This method creates the stages checkpoint list, either from a user-defined checkpoint list or,
        if a bool is provided, it creates a list with all the stages in the execution list.
        '''

        # Check if the checkpoint stages list has been already created
        if not hasattr(self, "__stages_to_checkpoint_list"):

            if self.settings["orchestrator"]["settings"]["stage_checkpoints"].IsBool():
                if self.settings["orchestrator"]["settings"]["stage_checkpoints"].GetBool():
                    # All stages are to be checkpointed
                    self.__stages_to_checkpoint_list = self.__GetExecutionList()
                else:
                    # Create an empty list to avoid errors
                    self.__stages_to_checkpoint_list = []

            elif self.settings["orchestrator"]["settings"]["stage_checkpoints"].IsStringArray():
                # Check that the stages asked to be checkpointed are in the execution list
                stages_asked_to_checkpoint = self.settings["orchestrator"]["settings"]["stage_checkpoints"].GetStringArray()
                aux_diff = set(stages_asked_to_checkpoint) - set(self.__GetExecutionList())
                if len(aux_diff) == 0:
                    self.__stages_to_checkpoint_list = stages_asked_to_checkpoint
                else:
                    err_msg = f"User asked to checkpoint stages that are not present in 'execution_list'. Please remove these {list(aux_diff)}."
                    raise ValueError(err_msg)
            else:
                err_msg = "'stage_checkpoints' value must be bool or a list with the stages to be checkpointed."
                raise TypeError(err_msg)
        
        return self.__stages_to_checkpoint_list
    
    def __GetStagesCheckpointFolder(self):
        '''Gets the folder to store the checkpoint files from the settings.'''

        if self.settings["orchestrator"]["settings"].Has("stage_checkpoints_folder"):
            return self.settings["orchestrator"]["settings"]["stage_checkpoints_folder"].GetString()
        else:
            warn_msg = "'stage_checkpoints_folder' is not provided. Creating a default 'checkpoints' one."
            KratosMultiphysics.Logger.PrintWarning(warn_msg)
            return 'checkpoints'
        
    def __OutputValidatedSettings(self):
        if KratosMultiphysics.ParallelEnvironment.GetDefaultDataCommunicator().Rank() == 0:
            output_validated_settings = False
            if self.settings["orchestrator"]["settings"].Has("output_validated_settings"):
                if self.settings["orchestrator"]["settings"]["output_validated_settings"].IsBool():
                    # If a bool is provided, we save the settings with a default name
                    output_validated_settings = self.settings["orchestrator"]["settings"]["output_validated_settings"].GetBool()
                    output_validated_settings_name = "ProjectParametersValidated.json"
                elif self.settings["orchestrator"]["settings"]["output_validated_settings"].IsString():
                    # If a string is provided it is taken as output name
                    output_validated_settings = True
                    output_validated_settings_name = self.settings["orchestrator"]["settings"]["output_validated_settings"].GetString()
                else:
                    err_msg = "'output_validated_settings' type is not supported. It is expected to be either a bool or a string."
                    raise TypeError(err_msg)

            if output_validated_settings:
                with open(f"{output_validated_settings_name}", 'w') as parameter_output_file:
                    parameter_output_file.write(self.settings.PrettyPrintJsonString())
