import os
import sys
import pickle
import pathlib
import importlib

import KratosMultiphysics

class Project:

    def __init__(self) -> None:
        '''
        
        Explain the members
        '''

        self.output_data = {}
        self.active_stages = {}
        self.model = KratosMultiphysics.Model()

    def GetModel(self):
        return self.model
    
    def AddActiveStage(self, stage_name, stage_instance):
        '''Adds the provided stage instance to the active stages dictionary.'''

        if self.active_stages.has_key(stage_name):
            err_msg = f"Stage '{stage_name}' is already active and cannot be added to active stages."
            raise Exception(err_msg)
        self.active_stages[stage_name] = stage_instance

    def RemoveActiveStage(self, stage_name):
        '''Removes an active stage instance from the current stages dictionary.'''

        del self.active_stages[stage_name]

    def Save(self, save_folder_name, file_name):
        '''Saves the Project current status.'''

        # Set the list of modules (Kratos and non-Kratos) that have been added up to current save
        required_modules = list(sys.modules.keys())

        # Create save folder
        checkpoint_path = pathlib.Path(save_folder_name)
        KratosMultiphysics.FilesystemExtensions.MPISafeCreateDirectories(checkpoint_path)

        # Save current status
        with open(os.path.join(checkpoint_path, file_name), 'wb+') as checkpoint_file:
            # Serialize current model and stages
            serializer = KratosMultiphysics.StreamSerializer()
            serializer.Save("Model", self.model)
            stage_names_list = []
            stage_instances_list = []
            for stage_name, stage_instance in self.active_stages.items():
                stage_instance.Save(serializer) # Make stage instance pickable by serializing and deleting all Kratos objects
                stage_names_list.append(stage_name) # Append current stage name
                stage_instances_list.append(stage_instance) # Append current stage pickable instance
            
            pickle.dump({
                "serializer" : serializer,
                "output_data" : self.output_data,
                "stage_names" : stage_names_list,
                "stage_instances" : stage_instances_list,
                "required_modules" : required_modules
            }, checkpoint_file, protocol=2)
   
    def Load(self, loading_point):
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
            loaded_data["serializer"].Load("Model", self.model)

            # Load stages
            for stage_name, stage_instance in zip(loaded_data["stage_names"], loaded_data["stage_instances"]):
                stage_instance.Load(loaded_data["serializer"]) # Take stage instance and restore it to its status before serialization
                self.active_stages[stage_name] = stage_instance # Append current stage instance to the active stages list

            # Load output data dictionary
            self.output_data = loaded_data["output_data"]
 