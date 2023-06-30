import abc
import sys
import importlib

import KratosMultiphysics
from KratosMultiphysics.model_parameters_factory import KratosModelParametersFactory

class MultistageOrchestrator(abc.ABC):

    def __init__(self, settings) -> None:
        self.settings = settings
        self.project = KratosMultiphysics.Project()

    @abc.abstractmethod
    def Run(self):
        '''Main function that runs the complete multistage simulation.'''

        pass

    def CheckStageSettings(self, stage_name):
        '''Check the settings for the given stage name
        
        This methods performs the check of the given stage name settings.
        Note that the Check() method of the stage instance is intentionally not called in here.
        This is due to the fact that the Check() of the stage instance is to be called within its Initialize() method.
        '''

        if self.settings["stages"][stage_name].Has("stage_postprocess"):
            if self.settings["stages"][stage_name]["stage_postprocess"].Has("modelers"):
                err_msg = f"Found 'modelers' field in 'stage_postprocess' of stage {stage_name}."
                err_msg += " Place the 'modelers' section in the next stage 'stage_preprocess'."
                raise Exception(err_msg)

    def CreateStage(self, stage_name):
        '''This method creates a stage instance

        Given a stage name, this method creates the corresponding stage instance from the "analysis_stage" registry entry in the settings.
        The stage instance is to be created from the registry. If not registered, it tries to use "analysis_stage" as Python module.
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
                err_msg = f"Please provide the 'analysis_stage' registry entry with module name (not the 'All' one) in '{stage_name}'."
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

    def RunCurrentStagePreprocess(self, stage_name):
        '''This function executes the preprocess of current stage.

        Note that the stage preprocess involves the execution of modelers and operations.
        '''

        if self.settings["stages"][stage_name].Has("stage_preprocess"):
            if self.settings["stages"][stage_name]["stage_preprocess"].Has("modelers"):
                for modeler in self.__GetModelers(stage_name):
                    modeler.SetupGeometryModel()
                for modeler in self.__GetModelers(stage_name):
                    modeler.PrepareGeometryModel()
                for modeler in self.__GetModelers(stage_name):
                    modeler.SetupModelPart()

            if self.settings["stages"][stage_name]["stage_preprocess"].Has("operations"):
                for operation in self.__GetOperations("stage_preprocess", stage_name):
                    operation.Execute()

    def RunCurrentStagePostprocess(self, stage_name):
        '''This function executes the postprocessing of current stage.
        
        Note that the stage postprocess deliberately involves operations only.
        '''

        if self.settings["stages"][stage_name].Has("stage_postprocess"):
            if self.settings["stages"][stage_name]["stage_postprocess"].Has("operations"):
                for operation in self.__GetOperations("stage_postprocess", stage_name):
                    operation.Execute()

    def GetProject(self):
        '''Returns the project.'''

        return self.project

    def GetNumberOfStages(self):
        '''Returns the number of stages.'''

        return len(self.GetExecutionList())

    def GetExecutionList(self):
        '''Creates and returns the execution list.
        This method creates the execution list, either from a user-defined execution list or,
        if this is not provided, from the stages declaration order.
        '''

        # Check if the execution list has been already created
        if not hasattr(self, "__execution_list"):
            # Default case in which the execution list is provided by the user
            if self.settings["project_settings"].Has("execution_list"):
                self.__execution_list = self.settings["project_settings"]["execution_list"].GetStringArray()
            # If not provided, create an auxiliary execution list from the stages declaration order
            else:
                KratosMultiphysics.Logger.PrintInfo("'execution_list' is not provided. Stages will be executed according to their declaration order.")
                self.__execution_list = list(self.settings["stages"].keys())

        return self.__execution_list
    
    def __GetModelers(self, stage_name):
        '''This method creates the modelers at the preprocess execution point.'''

        execution_point_settings = self.settings["stages"][stage_name]["stage_preprocess"]

        factory = KratosModelParametersFactory(self.model)
        return factory.ConstructListOfItems(execution_point_settings["modelers"])

    def __GetOperations(self, execution_point, stage_name):
        '''This method creates the operations at any execution point.'''

        if execution_point not in ["stage_preprocess","stage_postprocess"]:
            err_msg = f"Wrong execution point '{execution_point}'. Supported ones are 'stage_preprocess' and 'stage_postprocess'."
            raise Exception(err_msg)
        execution_point_settings = self.settings["stages"][stage_name][execution_point]

        factory = KratosModelParametersFactory(self.model)
        return factory.ConstructListOfItems(execution_point_settings["operations"])
