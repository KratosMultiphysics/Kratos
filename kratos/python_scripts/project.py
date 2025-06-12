import sys
import pickle
import pathlib
import importlib
from typing import Optional

import KratosMultiphysics
from KratosMultiphysics.analysis_stage import AnalysisStage

class Project:
    '''Kratos Multiphysics multistage project container.

    This class has two main purposes. First one is to hold the multistage components (output data, active stages and model)
    Second one is to perform the checkpoint save and load operations.

    Member variables:
    __settings -- Kratos parameters object with the multistage simulation settings
    __output_data -- Dictionary containing the stages data retrieved from GetFinalData
    __active_stages -- Dictionary containing the active (alive) stage instances
    __model -- Model instance
    '''

    def __init__(self, settings: KratosMultiphysics.Parameters) -> None:
        '''Constructs the multistage project container instance and sets current Kratos version in the settings.'''

        # Declare member variables
        self.__output_data: dict = {}
        self.__active_stages: dict = {}
        self.__settings: KratosMultiphysics.Parameters = settings
        self.__model: KratosMultiphysics.Model = KratosMultiphysics.Model()

        # Add Kratos version and compilation to settings
        kratos_version = f"{KratosMultiphysics.KratosGlobals.Kernel.Version()}-{KratosMultiphysics.KratosGlobals.Kernel.BuildType()}"
        self.__settings.AddString("kratos_version", kratos_version)

    def GetModel(self) -> KratosMultiphysics.Model:
        '''Returns the current multistage simulation model.'''

        return self.__model

    def GetSettings(self) -> KratosMultiphysics.Parameters:
        '''Returns the current multistage simulation settings.'''

        return self.__settings

    def GetOutputData(self) -> dict:
        '''Returns the current multistage simulation output data container.'''

        return self.__output_data

    def GetActiveStages(self) -> dict:
        '''Returns the current multistage simulation active stages dictionary.'''

        return self.__active_stages

    def AddActiveStage(self, stage_name: str, stage_instance: AnalysisStage) -> None:
        '''Adds the provided stage instance to the active stages dictionary.'''

        if self.__active_stages.has_key(stage_name):
            err_msg = f"Stage '{stage_name}' is already active and cannot be added to active stages."
            raise Exception(err_msg)
        self.__active_stages[stage_name] = stage_instance

    def RemoveActiveStage(self, stage_name : str) -> None:
        '''Removes an active stage instance from the current stages dictionary.'''

        del self.__active_stages[stage_name]

    def Save(self, save_folder_path: pathlib.Path, checkpoint_file_path: pathlib.Path, output_settings_file_path: Optional[pathlib.Path] = None) -> None:
        '''Saves the Project current status.'''

        # Set the list of modules (Kratos and non-Kratos) that have been added up to current save
        required_modules = list(sys.modules.keys())

        # Create save folder
        KratosMultiphysics.FilesystemExtensions.MPISafeCreateDirectories(save_folder_path)

        # Save current status
        with open(save_folder_path / checkpoint_file_path, 'wb+') as checkpoint_file:
            # Serialize current model and stages
            serializer = KratosMultiphysics.StreamSerializer()
            serializer_flags = KratosMultiphysics.Serializer.SHALLOW_GLOBAL_POINTERS_SERIALIZATION
            serializer.Set(serializer_flags)

            serializer.Save("SerializerFlags", serializer_flags) # Save the serializer flags (will be required for the loading)
            serializer.Save("Model", self.__model)
            stage_names_list = []
            stage_instances_list = []
            for stage_name, stage_instance in self.__active_stages.items():
                stage_instance.Save(serializer) # Make stage instance pickable by serializing and deleting all Kratos objects
                stage_names_list.append(stage_name) # Append current stage name
                stage_instances_list.append(stage_instance) # Append current stage pickable instance

            pickle.dump({
                "serializer" : serializer,
                "output_data" : self.__output_data,
                "stage_names" : stage_names_list,
                "stage_instances" : stage_instances_list,
                "required_modules" : required_modules
            }, checkpoint_file, protocol=2)

        # Output current settings
        if KratosMultiphysics.ParallelEnvironment.GetDefaultDataCommunicator().Rank() == 0:
            if output_settings_file_path:
                with open(save_folder_path / output_settings_file_path, 'w') as parameter_output_file:
                    parameter_output_file.write(self.__settings.PrettyPrintJsonString())

    def Load(self, loading_path: pathlib.Path) -> None:
        '''Loads a saved Project status into current one.'''

        # Load save path file
        with open(loading_path, 'rb') as loading_file:
            loaded_data = pickle.load(loading_file)

            # Load required modules
            # Note that this requires to be done first in order to have the variables in the model
            for module_name in loaded_data["required_modules"]:
                importlib.import_module(module_name)

            # Load and prepare serializer
            # This means to set the serializer flags before using it
            serializer = loaded_data["serializer"]
            serializer_flags = KratosMultiphysics.Flags()
            serializer.Load("SerializerFlags", serializer_flags)
            serializer.Set(serializer_flags)

            # Load model
            serializer.Load("Model", self.__model)

            # Load stages
            for stage_name, stage_instance in zip(loaded_data["stage_names"], loaded_data["stage_instances"]):
                stage_instance.Load(serializer) # Take stage instance and restore it to its status before serialization
                self.__active_stages[stage_name] = stage_instance # Append current stage instance to the active stages list

            # Load output data dictionary
            self.__output_data = loaded_data["output_data"]
