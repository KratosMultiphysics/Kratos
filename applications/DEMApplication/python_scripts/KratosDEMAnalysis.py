import KratosMultiphysics
import KratosMultiphysics.DEMApplication
from KratosMultiphysics.DEMApplication.DEM_analysis_stage import DEMAnalysisStage

import time
import sys
""" model = KratosMultiphysics.Model()
solution = Main.Solution(model)
solution.Run() """

class DEMAnalysisStageWithFlush(DEMAnalysisStage):

    def __init__(self, model, project_parameters, flush_frequency=10.0):
        super().__init__(model, project_parameters)
        self.flush_frequency = flush_frequency
        self.last_flush = time.time()

    def FinalizeSolutionStep(self):
        super().FinalizeSolutionStep()

        if self.parallel_type == "OpenMP":
            now = time.time()
            if now - self.last_flush > self.flush_frequency:
                sys.stdout.flush()
                self.last_flush = now

if __name__ == "__main__":
    from KratosMultiphysics import Logger
    Logger.GetDefaultOutput().SetSeverity(Logger.Severity.INFO)
    with open("ProjectParametersDEM.json",'r') as parameter_file:
        project_parameters = KratosMultiphysics.Parameters(parameter_file.read())

    model = KratosMultiphysics.Model()
    DEMAnalysisStageWithFlush(model, project_parameters).Run()
