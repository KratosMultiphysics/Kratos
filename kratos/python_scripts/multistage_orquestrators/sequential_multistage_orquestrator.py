import KratosMultiphysics
from KratosMultiphysics.multistage_orquestrators.multistage_orquestrator import MultistageOrquestrator

class SequentialMultistageOrquestrator(MultistageOrquestrator):

    def __init__(self, settings) -> None:
        super.__init__(settings)

    def Run(self):
        '''Main function that runs the complete multistage simulation.'''

        # Set the stages to be saved first
        # We do it first to warn the user about the non-present stages before running the first one
        self.__GetStagesToCheckpointList()

        # Check if loading from checkpoint is required and load
        if self.settings["project_settings"].Has("load_from_checkpoint"):
            loading_point = self.settings["project_settings"]["load_from_checkpoint"].GetString()
            self.GetProject().Load(loading_point)
            # self.__LoadCurrentStageCheckpoint()

        # Run the stages list
        for stage_name in self.GetExecutionList():
            # Check current stage settings
            self.CheckStageSettings(stage_name)

            # Set up and check current stage instance
            current_stage = self.CreateStage(stage_name)
            self.GetProject().AddActiveStage(stage_name, current_stage)

            # Execute current stage preprocess
            self.RunCurrentStagePreprocess(stage_name)

            # Run (solve) current stage
            current_stage.Run()

            # Execute current stage postprocess
            # Note that this is deliberately done before getting the final data in case this needs to be modified
            self.RunCurrentStagePostprocess(stage_name)

            # Get the final data dictionary
            self.GetProject().output_data[stage_name] = current_stage.GetAnalysisStageFinalData()

            # Check the current stage is to be checkpointed and save it if so
            if stage_name in self.__GetStagesToCheckpointList():
                save_folder_name = self.__GetStagesCheckpointFolder()
                self.GetProject().Save(save_folder_name, stage_name)
                # self.__SaveCurrentStageCheckpoint()
                
            # Output complete simulation settings (with defaults) as simulation report
            if (s:=self.settings["stages"][stage_name]) and s.Has("print_validated_settings") and s["print_validated_settings"].GetBool():
                with open(f"{stage_name}ValidatedSettings.json", 'w') as parameter_output_file:
                    parameter_output_file.write(current_stage.project_parameters.PrettyPrintJsonString())

            # Delete current stage instance
            self.GetProject().RemoveActiveStage(stage_name)

    def __GetStagesToCheckpointList(self):
        '''Creates and returns the stages checkpoint list.

        This method creates the stages checkpoint list, either from a user-defined checkpoint list or,
        if a bool is provided, it creates a list with all the stages in the execution list.
        '''

        # Check if the checkpoint stages list has been already created
        if not hasattr(self, "__stages_to_checkpoint_list"):

            if self.settings["project_settings"]["stage_checkpoints"].IsBool():
                if self.settings["project_settings"]["stage_checkpoints"].GetBool():
                    # All stages are to be checkpointed
                    self.__stages_to_checkpoint_list = self.GetExecutionList()
                else:
                    # Create an empty list to avoid errors
                    self.__stages_to_checkpoint_list = []

            elif self.settings["project_settings"]["stage_checkpoints"].IsStringArray():
                # Check that the stages asked to be checkpointed are in the execution list
                stages_asked_to_checkpoint = self.settings["project_settings"]["stage_checkpoints"].GetStringArray()
                aux_diff = set(stages_asked_to_checkpoint) - set(self.GetExecutionList())
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

        if self.settings["project_settings"].Has("stage_checkpoints_folder"):
            return self.settings["project_settings"]["stage_checkpoints_folder"].GetString()
        else:
            warn_msg = "'stage_checkpoints_folder' is not provided. Creating a default 'checkpoints' one."
            KratosMultiphysics.Logger.PrintWarning(warn_msg)
            return 'checkpoints'
