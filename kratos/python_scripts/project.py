import os
import sys
import pickle
import pathlib
import importlib

import KratosMultiphysics

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

    def __init__(self, settings : KratosMultiphysics.Parameters) -> None:
        '''Constructs the multistage project container instance and sets current Kratos version in the settings.'''

        # Declare member variables
        self.__output_data = {}
        self.__active_stages = {}
        self.__settings = settings
        self.__model = KratosMultiphysics.Model()

        # Add Kratos version and compilation to settings
        kratos_version = f"{KratosMultiphysics.KratosGlobals.Kernel.Version()}-{KratosMultiphysics.KratosGlobals.Kernel.BuildType()}"
        self.__settings.AddString("kratos_version", kratos_version)

    def GetModel(self):
        '''Returns the current multistage simulation model.'''

        return self.__model
    
    def GetSettings(self):
        '''Returns the current multistage simulation settings.'''

        return self.__settings
    
    def GetOutputData(self):
        '''Returns the current multistage simulation output data container.'''

        return self.__output_data
    
    def GetActiveStages(self):
        '''Returns the current multistage simulation active stages dictionary.'''

        return self.__active_stages
    
    def AddActiveStage(self, stage_name : str, stage_instance):
        '''Adds the provided stage instance to the active stages dictionary.'''

        if self.__active_stages.has_key(stage_name):
            err_msg = f"Stage '{stage_name}' is already active and cannot be added to active stages."
            raise Exception(err_msg)
        self.__active_stages[stage_name] = stage_instance

    def RemoveActiveStage(self, stage_name : str):
        '''Removes an active stage instance from the current stages dictionary.'''

        del self.__active_stages[stage_name]

    def Save(self, save_folder_name : str, checkpoint_file_name : str, output_settings_file_name : str = None):
        '''Saves the Project current status.'''

        # Set the list of modules (Kratos and non-Kratos) that have been added up to current save
        required_modules = list(sys.modules.keys())

        # Create save folder
        checkpoint_path = pathlib.Path(save_folder_name)
        KratosMultiphysics.FilesystemExtensions.MPISafeCreateDirectories(checkpoint_path)

        # Save current status
        with open(os.path.join(checkpoint_path, checkpoint_file_name), 'wb+') as checkpoint_file:
            # Serialize current model and stages
            serializer = KratosMultiphysics.StreamSerializer()
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
            if output_settings_file_name:
                with open(os.path.join(checkpoint_path, output_settings_file_name), 'w') as parameter_output_file:
                    parameter_output_file.write(self.__settings.PrettyPrintJsonString())
   
    def Load(self, loading_point : str):
        '''Loads a saved Project status into current one.'''

        # Load save path file
        loading_path = pathlib.Path(loading_point)
        with open(loading_path, 'rb') as loading_file:
            loaded_data = pickle.load(loading_file)

            # Load required modules
            # Note that this requires to be done first in order to have the variables in the model
            for module_name in loaded_data["required_modules"]:
                importlib.import_module(module_name)

            # Load model
            loaded_data["serializer"].Load("Model", self.__model)

            # Load stages
            for stage_name, stage_instance in zip(loaded_data["stage_names"], loaded_data["stage_instances"]):
                stage_instance.Load(loaded_data["serializer"]) # Take stage instance and restore it to its status before serialization
                self.__active_stages[stage_name] = stage_instance # Append current stage instance to the active stages list

            # Load output data dictionary
            self.__output_data = loaded_data["output_data"]
 