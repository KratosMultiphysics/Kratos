import sys
import time
import importlib

import KratosMultiphysics
from KratosMultiphysics.modeler_factory import KratosModelerFactory
from KratosMultiphysics.model_parameters_factory import KratosModelParametersFactory

class MultistageAnalysis(object):

    def __init__(self, model, project_parameters) -> None:
        self.model = model
        self.current_stage_name = None
        self.settings = project_parameters
        self.stages_map = self.__CreateStagesMap()

    def Check(self):
        for stage_name in self.settings["execution_list"].GetStringArray():
            self.CheckStage(stage_name)

    def CheckStage(self, stage_name):
        if self.settings["stages"][stage_name].Has("stage_postprocess"):
            if self.settings["stages"][stage_name]["stage_postprocess"].Has("modelers"):
                err_msg = f"Found 'modelers' field in 'stage_postprocess' of stage {stage_name}."
                err_msg += f" Place the 'modelers' section in the next stage 'stage_postprocess'."
                raise Exception(err_msg)

    def Run(self):
        # First check all the stages input
        self.Check()

        # Run the stages list
        for stage_name in self.settings["execution_list"].GetStringArray():
            self.current_stage_name = stage_name
            self.RunCurrentStagePreprocess()
            self.RunCurrentStage()
            self.RunCurrentStagePostprocess()

    def RunCurrentStagePreprocess(self):
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
        self.GetCurrentStage().Run()

    def RunCurrentStagePostprocess(self):
        if self.settings["stages"][self.GetCurrentStageName()].Has("stage_postprocess"):
            if self.settings["stages"][self.GetCurrentStageName()]["stage_postprocess"].Has("operations"):
                for operation in self.__GetOperations("stage_postprocess"):
                    operation.Execute()

    def GetNumberOfStages(self):
        return len(self.stages_map)

    def GetCurrentStage(self):
        return self.stages_map[self.current_stage_name]

    def GetCurrentStageName(self):
        return self.current_stage_name

    def __GetModelers(self):
        #TODO: Add error thrown for wrong execution_point argument
        execution_point_settings = self.settings["stages"][self.GetCurrentStageName()]["stage_preprocess"]

        factory = KratosModelParametersFactory(self.model)
        return factory.ConstructListOfItems(execution_point_settings["modelers"])

    def __GetOperations(self, execution_point):
        #TODO: Add error thrown for wrong execution_point argument
        execution_point_settings = self.settings["stages"][self.GetCurrentStageName()][execution_point]

        factory = KratosModelParametersFactory(self.model)
        return factory.ConstructListOfItems(execution_point_settings["operations"])

    def __CreateStagesMap(self):
        stages_map = {}
        for stage_name in self.settings["execution_list"].GetStringArray():
            analysis_stage_module_name = self.settings["stages"][stage_name]["analysis_stage"].GetString()
            analysis_stage_class_name = analysis_stage_module_name.split('.')[-1]
            analysis_stage_class_name = ''.join(x.title() for x in analysis_stage_class_name.split('_'))

            analysis_stage_module = importlib.import_module(analysis_stage_module_name)
            analysis_stage_class = getattr(analysis_stage_module, analysis_stage_class_name)

            stages_map[stage_name] = self.__CreateAnalysisStageWithFlushInstance(analysis_stage_class, self.model,  KratosMultiphysics.Parameters(self.settings["stages"][stage_name]))
        return stages_map

    #TODO: Think about this
    def __CreateAnalysisStageWithFlushInstance(self, cls, global_model, parameters):
        class AnalysisStageWithFlush(cls):

            def __init__(self, model,project_parameters, flush_frequency=10.0):
                super().__init__(model,project_parameters)
                self.flush_frequency = flush_frequency
                self.last_flush = time.time()
                sys.stdout.flush()

            def Initialize(self):
                super().Initialize()
                sys.stdout.flush()

            def FinalizeSolutionStep(self):
                super().FinalizeSolutionStep()

                if self.parallel_type == "OpenMP":
                    now = time.time()
                    if now - self.last_flush > self.flush_frequency:
                        sys.stdout.flush()
                        self.last_flush = now

            def __repr__(self) -> str:
                return super().__repr__()

        return AnalysisStageWithFlush(global_model, parameters)
