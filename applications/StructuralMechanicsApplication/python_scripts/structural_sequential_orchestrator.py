import pathlib
from typing import Dict, Optional
import KratosMultiphysics
from KratosMultiphysics.project import Project
from KratosMultiphysics.orchestrators.sequential_orchestrator import SequentialOrchestrator

class StructuralSequentialOrchestrator(SequentialOrchestrator):

    def CreateStage(self, stage_name: str):
        self.current_stage = super().CreateStage(stage_name)
        return self.current_stage
    
    def RunCurrentStagePreprocess(self, stage_name: str, data: Optional[Dict] = None):
        super().RunCurrentStagePreprocess(stage_name, data)
        if self.GetProject().GetSettings()["stages"][stage_name]["stage_preprocess"].Has("get_output_data"):
            output_stage_name = self.GetProject().GetSettings()["stages"][stage_name]["stage_preprocess"]["get_output_data"].GetString()
            self.current_stage.GetProjectOutputData(output_stage_name, self.GetProject().GetOutputData())