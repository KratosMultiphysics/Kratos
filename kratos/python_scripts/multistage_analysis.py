import importlib

import KratosMultiphysics
from KratosMultiphysics.model_parameters_factory import KratosModelParametersFactory

class MultistageAnalysis:

    def __init__(self, model, project_parameters) -> None:
        self.model = model
        self.__current_stage_name = None
        self.settings = project_parameters
        self.__stages_map = self.__CreateStagesMap()

    def Check(self):
        '''Performs the check of the complete multistage simulation.'''

        for stage_name in self.__GetExecutionList():
            self.CheckStage(stage_name)

    def CheckStage(self, stage_name):
        '''Performs the check of a single stage from a multistage simulation.'''

        if self.settings["stages"][stage_name].Has("stage_postprocess"):
            if self.settings["stages"][stage_name]["stage_postprocess"].Has("modelers"):
                err_msg = f"Found 'modelers' field in 'stage_postprocess' of stage {stage_name}."
                err_msg += " Place the 'modelers' section in the next stage 'stage_preprocess'."
                raise Exception(err_msg)

    def Run(self):
        '''Main function that runs the complete multistage simulation.'''

        # First check all the stages input
        self.Check()

        # Run the stages list
        for stage_name in self.__GetExecutionList():
            self.__current_stage_name = stage_name
            self.RunCurrentStagePreprocess()
            self.RunCurrentStage()
            self.RunCurrentStagePostprocess()

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

    def RunCurrentStage(self):
        '''This function executes (solves) the current stage.
        Note that this call is equivalent to the traditional single-stage simulation run call.
        '''

        self.GetCurrentStage().Run()

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

        return self.__stages_map[self.__current_stage_name]

    def GetCurrentStageName(self):
        '''Returns the current stage name.'''
        return self.__current_stage_name

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

    def __CreateStagesMap(self):
        '''This method creates the stages map.
        The keys of the map are the stage names, which are retrieved from the execution list (see __GetExecutionList()).
        The values of the map are the stages instances, which are created by importing its corresponding analysis stage Python module.
        '''

        stages_map = {}
        for stage_name in self.__GetExecutionList():
            # Get the Python module name and class name implementing each analysis stage
            # Note that we assume that the class name is the provided module (file) name in CamelCase
            analysis_stage_module_name = self.settings["stages"][stage_name]["analysis_stage"].GetString()
            analysis_stage_class_name = analysis_stage_module_name.split('.')[-1]
            analysis_stage_class_name = ''.join(x.title() for x in analysis_stage_class_name.split('_'))

            # Import the stage module and create the corresponding instance
            analysis_stage_module = importlib.import_module(analysis_stage_module_name)
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

            # Insert current stage instance in the stages map
            stages_map[stage_name] = stage_instance

        return stages_map
