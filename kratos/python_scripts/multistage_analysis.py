import os
import sys
import pickle
import pathlib
import importlib

import KratosMultiphysics
from KratosMultiphysics.model_parameters_factory import KratosModelParametersFactory

class MultistageAnalysis:

    def __init__(self, project_parameters) -> None:
        self.output_data = {}
        self.settings = project_parameters
        self.model = KratosMultiphysics.Model()
        self.initial_modules = [k for k in sys.modules.keys()]

    def CheckStage(self, stage_name):
        '''Performs the check of a single stage from a multistage simulation.'''

        if self.settings["stages"][stage_name].Has("stage_postprocess"):
            if self.settings["stages"][stage_name]["stage_postprocess"].Has("modelers"):
                err_msg = f"Found 'modelers' field in 'stage_postprocess' of stage {stage_name}."
                err_msg += " Place the 'modelers' section in the next stage 'stage_preprocess'."
                raise Exception(err_msg)

    def Run(self):
        '''Main function that runs the complete multistage simulation.'''

        # Check if loading from checkpoint is required and load
        if self.settings["checkpoint_settings"]["load"].GetBool():
            self.__LoadCurrentStageCheckpoint()

        # Run the stages list
        for stage_name in self.__GetExecutionList():
            # Set up and check current stage instance
            self.SetCurrentStageName(stage_name)
            current_stage = self.GetCurrentStage()
            self.CheckStage(self.GetCurrentStageName())

            # Execute current stage preprocess
            self.RunCurrentStagePreprocess()

            # Run (solve) current stage
            current_stage.Run()

            # Execute current stage postprocess
            # Note that this is deliberately done before getting the final data in case this needs to be modified
            self.RunCurrentStagePostprocess()

            # Get the final data dictionary
            self.output_data[self.GetCurrentStageName()] = current_stage.GetFinalData() #NOTE: only python types in here - TODO: can we check that?

            # Check the current stage is to be checkpointed and save it if so
            if self.GetCurrentStageName() in self.__GetStagesToCheckpointList():
                self.__SaveCurrentStageCheckpoint()
                
            # self.print_project_paramters_report(self.GetCurrentStageName()) //FIXME: hay un prettyprintjsonstring

            # Delete current stage instance
            del current_stage

    def __SaveCurrentStageCheckpoint(self):
        # Set the list of Kratos modules that have been added up to current stage
        required_modules = [m for m in (sys.modules.keys() - self.initial_modules) if m.split(".")[0] == "KratosMultiphysics"]

        # Create checkpoints folder
        checkpoint_path = pathlib.Path(self.settings["checkpoint_settings"]["save_folder_name"].GetString())
        KratosMultiphysics.FilesystemExtensions.MPISafeCreateDirectories(checkpoint_path)

        # Save current checkpoint
        with open(os.path.join(checkpoint_path, self.GetCurrentStageName()), 'wb+') as checkpoint_file:
            serializer = KratosMultiphysics.StreamSerializer()
            serializer.Save("Model",self.model)
            pickle.dump({
                "model_serializer" : serializer,
                "output_data" : self.output_data,
                "required_modules" : required_modules
            }, checkpoint_file, 2)

    def __LoadCurrentStageCheckpoint(self):
        checkpoint_path = pathlib.Path(self.settings["checkpoint_settings"]["load_checkpoint"].GetString())
        with open(checkpoint_path, 'rb') as checkpoint_file:
            loaded_data = pickle.load(checkpoint_file)

            # Load required modules
            # Note that this requires to be done first in order to have the variables in the model
            for module_name in loaded_data["required_modules"]:
                importlib.import_module(module_name)

            # Load model
            loaded_data["model_serializer"].Load("Model",self.model)

            # Load output data dictionary
            self.output_data = loaded_data["output_data"]

    def RunCurrentStagePreprocess(self):
        '''This function executes the preprocess of current stage.
        Note that the stage preprocess involves the execution of modelers and operations.
        '''

        if self.settings["stages"][self.GetCurrentStageName()].Has("stage_preprocess"):
            if self.settings["stages"][self.GetCurrentStageName()]["stage_preprocess"].Has("modelers"):
                for modeler in self.__GetModelers():
                    modeler.SetupGeometryModel()
                for modeler in self.__GetModelers():
                    modeler.PrepareGeometryModel()
                for modeler in self.__GetModelers():
                    modeler.SetupModelPart()

            if self.settings["stages"][self.GetCurrentStageName()]["stage_preprocess"].Has("operations"):
                for operation in self.__GetOperations("stage_preprocess"):
                    operation.Execute()

    def RunCurrentStagePostprocess(self):
        '''This function executes the postprocessing of current stage.
        Note that the stage postprocess deliberately involves operations only.
        '''

        if self.settings["stages"][self.GetCurrentStageName()].Has("stage_postprocess"):
            if self.settings["stages"][self.GetCurrentStageName()]["stage_postprocess"].Has("operations"):
                for operation in self.__GetOperations("stage_postprocess"):
                    operation.Execute()

    def GetNumberOfStages(self):
        '''Returns the number of stages.'''

        return len(self.__stages_map)

    def GetCurrentStage(self):
        '''Returns the current stage instance.'''

        return self.__CreateStage(self.GetCurrentStageName())

    def GetCurrentStageName(self):
        '''Gets the current stage name.'''

        return self.__current_stage_name

    def SetCurrentStageName(self, stage_name):
        '''Sets the current stage name.'''

        self.__current_stage_name = stage_name

    def __GetExecutionList(self):
        '''Creates and returns the execution list.
        This method creates the execution list, either from a user-defined execution list or,
        if this is not provided, from the stages declaration order.
        '''

        # Check if the execution list has been already created
        if not hasattr(self, "__execution_list"):
            # Default case in which the execution list is provided by the user
            if self.settings.Has("execution_list"):
                self.__execution_list = self.settings["execution_list"].GetStringArray()
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

            if self.settings["checkpoint_settings"]["save"].IsBool():
                if self.settings["checkpoint_settings"]["save"].GetBool():
                    # All stages are to be checkpointed
                    self.__stages_to_checkpoint_list = self.__GetExecutionList()
                else:
                    # Create an empty list to avoid errors
                    self.__stages_to_checkpoint_list = []

            elif self.settings["checkpoint_settings"]["save"].IsStringArray():
                # Check that the stages asked to be checkpointed are in the execution list
                stages_asked_to_checkpoint = self.settings["checkpoint_settings"]["save"].GetStringArray()
                aux_diff = stages_asked_to_checkpoint - self.__GetExecutionList()
                if len(aux_diff) == 0:
                    self.__stages_to_checkpoint_list = stages_asked_to_checkpoint
                else:
                    err_msg = f"User asked to checkpoint stages that are not present in 'execution_list'. Please remove these {list(aux_diff):?}."
                    raise ValueError(err_msg)
            else:
                err_msg = "'save' value must be bool or a list with the stages to be checkpointed."
                raise TypeError(err_msg)
        
        return self.__stages_to_checkpoint_list

    def __GetModelers(self):
        '''This method creates the modelers at the preprocess execution point.'''

        execution_point_settings = self.settings["stages"][self.GetCurrentStageName()]["stage_preprocess"]

        factory = KratosModelParametersFactory(self.model)
        return factory.ConstructListOfItems(execution_point_settings["modelers"])

    def __GetOperations(self, execution_point):
        '''This method creates the operations at any execution point.'''

        if execution_point not in ["stage_preprocess","stage_postprocess"]:
            err_msg = f"Wrong execution point '{execution_point}'. Supported ones are 'stage_preprocess' and 'stage_postprocess'."
            raise Exception(err_msg)
        execution_point_settings = self.settings["stages"][self.GetCurrentStageName()][execution_point]

        factory = KratosModelParametersFactory(self.model)
        return factory.ConstructListOfItems(execution_point_settings["operations"])

    def __CreateStage(self, stage_name):
        '''This method creates the stages map.
        The keys of the map are the stage names, which are retrieved from the execution list (see __GetExecutionList()).
        The values of the map are the stages instances, which are created by importing its corresponding analysis stage Python module.
        '''

        # Get the current analysis stage class and module names
        input_analysis_stage = self.settings["stages"][stage_name]["analysis_stage"].GetString()

        # Check if the current stage is to be obtained from the Kratos registry
        if input_analysis_stage.split(".")[0] == "Stages":
            # Import current stage main module to populate the registry
            # This is required as the applications Python stuff is registered when the application is imported first
            # Note that in here we are following the convention that the middle strings are the module name
            main_module_name = '.'.join(input_analysis_stage.split(".")[1:-1])
            if main_module_name == "All":
                err_msg = "Please provide the registry entry with module name (not the 'All' one)."
                raise ValueError(err_msg)
            try:
                importlib.import_module(main_module_name)
            except ImportError:
                err_msg = f"Unbable to import '{main_module_name}'."

            # Once we have populated the registry we retrieve the analysis stage from it
            # Note that in here we are following the registry convention in @python_registry_utilities.py
            if KratosMultiphysics.Registry.HasItem(input_analysis_stage):
                # Get current analysis data from registry
                analysis_stage_reg_entry = KratosMultiphysics.Registry[input_analysis_stage]
                analysis_stage_class_name = analysis_stage_reg_entry["ClassName"]
                if "ModuleName" in analysis_stage_reg_entry:
                    # Standard case in which the analysis stage is registered in a module
                    analysis_stage_module_name = analysis_stage_reg_entry["ModuleName"]
                else:
                    # Alternative in which only the "ClassName" is provided with no module (e.g. custom user-defined stages in MainKratos.py)
                    analysis_stage_module_name = None
            else:
                err_msg = f"Trying to retrieve a non-registered analysis stage '{input_analysis_stage}'"
                raise ValueError(err_msg)
        else:
            # Throw a warning to let the user know that the current analysis stage is not registered
            warn_msg = f"Current analysis stage '{input_analysis_stage}' is not registered. Assuming that provided string contains module and class names."
            KratosMultiphysics.Logger.PrintWarning(warn_msg)

            # Get the Python module name and class name implementing each analysis stage
            # Note that we assume that the class name is the provided module (file) name in CamelCase
            analysis_stage_module_name = input_analysis_stage
            analysis_stage_class_name = analysis_stage_module_name.split('.')[-1]
            analysis_stage_class_name = ''.join(x.title() for x in analysis_stage_class_name.split('_'))

        # Import the stage module
        # Note that for the case of custom analysis stage (i.e. no Kratos module) we check if the class exists in the __main__ module
        if analysis_stage_module_name is not None:
            analysis_stage_module = importlib.import_module(analysis_stage_module_name)
        else:
            analysis_stage_module = sys.modules["__main__"]
        
        # Create the analysis stage instance
        if hasattr(analysis_stage_module, analysis_stage_class_name):
            # First we check for the expected class name
            analysis_stage_class = getattr(analysis_stage_module, analysis_stage_class_name)
            stage_instance = analysis_stage_class(self.model,  KratosMultiphysics.Parameters(self.settings["stages"][stage_name]))
        elif hasattr(analysis_stage_module, "Create"):
            # If Kratos convention is not fulfilled we search for a Create method
            stage_instance = analysis_stage_module.Create(self.model,  KratosMultiphysics.Parameters(self.settings["stages"][stage_name]))
        else:
            err_msg = f"Analysis stage in '{analysis_stage_module_name}' Python module cannot be created. Please check class name or provide a 'Create' method."
            raise Exception(err_msg)
        
        return stage_instance
