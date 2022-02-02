import sys
import time
import importlib

import KratosMultiphysics

class MultistageAnalysis(object):

    def __init__(self, model, project_parameters) -> None:
        self.model = model
        self.current_stage_index = None
        self.settings = project_parameters
        self.stages_list = self.__CreateStagesList()

    def Run(self):
        for i_stage in range(self.GetNumberOfStages()):
            self.current_stage_index = i_stage
            self.InitializeCurrentStage()
            self.RunCurrentStage()
            self.FinalizeCurrentStage()

    def InitializeCurrentStage(self):
        pass

    def RunCurrentStage(self):
        self.stages_list[self.current_stage_index].Run()

    def FinalizeCurrentStage(self):
        pass

    def GetNumberOfStages(self):
        return len(self.stages_list)

    def GetCurrentStage(self):
        return self.stages_list[self.current_stage_index]

    def GetCurrentStageIndex(self):
        return self.current_stage_index

    def __CreateStagesList(self):
        stages_list = []
        for stage_settings in self.settings["stages"]:
            analysis_stage_module_name = stage_settings["analysis_stage"].GetString()
            analysis_stage_class_name = analysis_stage_module_name.split('.')[-1]
            analysis_stage_class_name = ''.join(x.title() for x in analysis_stage_class_name.split('_'))

            analysis_stage_module = importlib.import_module(analysis_stage_module_name)
            analysis_stage_class = getattr(analysis_stage_module, analysis_stage_class_name)

            stages_list.append(self.__CreateAnalysisStageWithFlushInstance(analysis_stage_class, self.model,  KratosMultiphysics.Parameters(stage_settings)))
        return stages_list

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
