import pathlib

import KratosMultiphysics
from KratosMultiphysics.project import Project
from KratosMultiphysics.orchestrators.orchestrator import Orchestrator

class SequentialOrchestrator(Orchestrator):
    '''Multistage orchestrator to sequentially execute a collection of stages.

    This multistage orchestrator executes the different stages sequentially (i.e. one after the other).
    The execution order is specified by the 'execution_list' in the input settings.
    The save and load (aka checkpointing) is executed between stages.

    Member variables:
    _execution_list -- list containing the names of the stages to be executed
    _stages_to_checkpoint_list -- list containing the names of the stages to be checkpointed (saved)
    '''

    def __init__(self, project : Project) -> None:
        '''Constructs the sequential multistage orchestrator instance.'''

        super().__init__(project)

    def Run(self) -> None:
        '''Main function that runs a sequential multistage simulation.'''

        # Set the stages to be saved first
        # We do it first to warn the user about the non-present stages before running the first one
        self.__GetStagesToCheckpointList()

        # Check if loading from checkpoint is required and load
        if self.GetProject().GetSettings()["orchestrator"]["settings"].Has("load_from_checkpoint"):
            if not self.GetProject().GetSettings()["orchestrator"]["settings"]["load_from_checkpoint"].IsNull():
                if not self.GetProject().GetSettings()["orchestrator"]["settings"]["load_from_checkpoint"].IsString():
                    err_msg = "'load_from_checkpoint' is expected to be null or a string with the loading point."
                    raise TypeError(err_msg)
                else:
                    loading_point = self.GetProject().GetSettings()["orchestrator"]["settings"]["load_from_checkpoint"].GetString()
                    self.GetProject().Load(pathlib.Path(loading_point))

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
            self.GetProject().GetOutputData()[stage_name] = current_stage.GetFinalData()

            # Execute current stage postprocess
            self.RunCurrentStagePostprocess(stage_name)

            # Delete current stage instance
            del current_stage

            # Output validated simulation settings (with defaults) as simulation report
            output_settings_name = self.__PrepareOutputSettings(stage_name)

            # Check the current stage is to be checkpointed and save it if so
            if stage_name in self.__GetStagesToCheckpointList():
                save_folder_name = self.__GetStagesCheckpointFolder()
                self.GetProject().Save(pathlib.Path(save_folder_name), pathlib.Path(stage_name), output_settings_name)


    def GetNumberOfStages(self) -> int:
        '''Returns the number of stages.'''

        return len(self.__GetExecutionList())

    def __GetExecutionList(self) -> list:
        '''Creates and returns the execution list.

        The first time it is called, this method creates the execution list,
        either from a user-defined execution list or, if this is not provided,
        from the stages declaration order.
        The subsequent calls directly return the already existent '_execution_list'.
        '''

        # Check if the execution list has been already created
        if not hasattr(self, "_execution_list"):
            # Default case in which the execution list is provided by the user
            if self.GetProject().GetSettings()["orchestrator"]["settings"].Has("execution_list"):
                self._execution_list = self.GetProject().GetSettings()["orchestrator"]["settings"]["execution_list"].GetStringArray()
            # If not provided, create an auxiliary execution list from the stages declaration order
            else:
                KratosMultiphysics.Logger.PrintInfo("'execution_list' is not provided. Stages will be executed according to their declaration order.")
                self._execution_list = list(self.GetProject().GetSettings()["stages"].keys())

        return self._execution_list

    def __GetStagesToCheckpointList(self) -> list:
        '''Creates and returns the stages checkpoint list.

        The first time it is called, this method creates the stages checkpoint list,
        either from a user-defined checkpoint list or, if a bool is provided, it creates
        a list with all the stages in the execution list.
        The subsequent calls directly return the already existent '_stages_to_checkpoint_list'.
        '''

        # Check if the checkpoint stages list has been already created
        if not hasattr(self, "_stages_to_checkpoint_list"):
            if self.GetProject().GetSettings()["orchestrator"]["settings"]["stage_checkpoints"].IsBool():
                if self.GetProject().GetSettings()["orchestrator"]["settings"]["stage_checkpoints"].GetBool():
                    # All stages are to be checkpointed
                    self._stages_to_checkpoint_list = self.__GetExecutionList()
                else:
                    # Create an empty list to avoid errors
                    self._stages_to_checkpoint_list = []

            elif self.GetProject().GetSettings()["orchestrator"]["settings"]["stage_checkpoints"].IsStringArray():
                # Check that the stages asked to be checkpointed are in the execution list
                stages_asked_to_checkpoint = self.GetProject().GetSettings()["orchestrator"]["settings"]["stage_checkpoints"].GetStringArray()
                aux_diff = set(stages_asked_to_checkpoint) - set(self.__GetExecutionList())
                if not aux_diff:
                    self._stages_to_checkpoint_list = stages_asked_to_checkpoint
                else:
                    err_msg = f"User asked to checkpoint stages that are not present in 'execution_list'. Please remove these {list(aux_diff)}."
                    raise ValueError(err_msg)
            else:
                err_msg = "'stage_checkpoints' value must be bool or a list with the stages to be checkpointed."
                raise TypeError(err_msg)

        return self._stages_to_checkpoint_list

    def __GetStagesCheckpointFolder(self) -> str:
        '''Gets the folder to store the checkpoint files from the settings.'''

        if self.GetProject().GetSettings()["orchestrator"]["settings"].Has("stage_checkpoints_folder"):
            return self.GetProject().GetSettings()["orchestrator"]["settings"]["stage_checkpoints_folder"].GetString()
        else:
            info_msg = "'stage_checkpoints_folder' is not provided. Creating a default 'checkpoints' one."
            KratosMultiphysics.Logger.PrintInfo("SequentialOrchestrator", info_msg)
            return 'checkpoints'

    def __PrepareOutputSettings(self, stage_name: str) -> str:
        '''Prepares the settings to be output.

        This function modifies current settings in order to prepare them to be output.
        This includes setting the loading checkpoint as the one to be saved and modifying the execution list accordingly.
        '''

        output_validated_settings = False
        output_validated_settings_name = None
        orchestrator_settings = self.GetProject().GetSettings()["orchestrator"]["settings"]
        if orchestrator_settings.Has("output_validated_settings"):
            if orchestrator_settings["output_validated_settings"].IsBool():
                # If a bool is provided, we save the settings with a default name
                output_validated_settings = orchestrator_settings["output_validated_settings"].GetBool()
                output_validated_settings_name = f"ProjectParametersValidated_{stage_name}.json"
            elif orchestrator_settings["output_validated_settings"].IsString():
                # If a string is provided it is taken as output name
                output_validated_settings = True
                user_defined_name = orchestrator_settings["output_validated_settings"].GetString()
                output_validated_settings_name = f"{user_defined_name}_{stage_name}.json"
            else:
                err_msg = "'output_validated_settings' type is not supported. It is expected to be either a bool or a string."
                raise TypeError(err_msg)

        # Setting as the next stage as default loading point
        if output_validated_settings:
            # Get next stage name
            # Note that we check if current stage is last one. If so, we keep current loading checkpoint value
            position = self.__GetExecutionList().index(stage_name)
            if not position == self.GetNumberOfStages() - 1:
                # Set the next stage as loading point
                loading_point = f"{self.__GetStagesCheckpointFolder()}/{stage_name}"
                if orchestrator_settings.Has("load_from_checkpoint"):
                    orchestrator_settings["load_from_checkpoint"].SetString(loading_point)
                else:
                    orchestrator_settings.AddString("load_from_checkpoint", loading_point)

                # Remove current checkpointed stage from the execution list
                # Note that we deliberately do it in a copy to avoid modifying the current execution one
                new_execution_list = self.__GetExecutionList().copy()
                new_execution_list.remove(stage_name)
                orchestrator_settings["execution_list"].SetStringArray(new_execution_list)
            else:
                # If there is load checkpoint remove it as this is the last stage
                if orchestrator_settings.Has("load_from_checkpoint"):
                    orchestrator_settings.RemoveValue("load_from_checkpoint")

                # Set an empty execution list
                new_execution_list = []
                orchestrator_settings["execution_list"].SetStringArray(new_execution_list)

        return output_validated_settings_name

