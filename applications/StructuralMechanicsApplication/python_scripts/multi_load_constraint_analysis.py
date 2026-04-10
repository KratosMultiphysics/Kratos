import KratosMultiphysics
from KratosMultiphysics.analysis_stage import AnalysisStage

class MultiLoadConstraintAnalysis(AnalysisStage):
    def __init__(self, model, project_parameters):
        self.model = model
        self.settings = project_parameters["solver_settings"]
        self.model_part_name = self.settings["model_part_name"].GetString()

        if self.model_part_name == "":
            raise Exception('Please specify a model_part name!')
        
        #check if model already has model part, otherwise create it and set domain size
        if self.model.HasModelPart(self.model_part_name):
            self.main_model_part = self.model[self.model_part_name]
        else:
            self.main_model_part = self.model.CreateModelPart(self.model_part_name)
            domain_size = self.settings["domain_size"].GetInt()
            if domain_size < 0:
                raise Exception('Please specify a "domain size" >= 0!')
            self.main_model_part.ProcessInfo.SetValue(KratosMultiphysics.DOMAIN_SIZE, domain_size)

        super().__init__(model, project_parameters)

    def Run(self):
        pass

    def GetProjectOutputData(self, stage_name, output_data):
        self.lhs = output_data[stage_name]["lhs"]
        self.rhss = output_data[stage_name]["rhss"]
        print(self.lhs)
        for key, val in self.rhss.items():
            print(key)
            print(val)