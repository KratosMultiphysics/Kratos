import KratosMultiphysics
import numpy as np
import KratosMultiphysics.MeshingApplication as KratosMA

def Factory(settings, model):
    if not isinstance(settings, KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return InterpolateVelocityProcess(model, settings["Parameters"])

## All the processes python should be derived from "Process"
class InterpolateVelocityProcess(KratosMultiphysics.Process):
    def __init__(self, model, settings):
        KratosMultiphysics.Process.__init__(self)

        default_settings = KratosMultiphysics.Parameters("""{
            "origin_model_part_file_name" : "NameOfMDPAfile",
            "destination_model_part_name" : "ModelPartName",
            "interpolation_variables_list": []
        }""")

        settings.ValidateAndAssignDefaults(default_settings)

        self.model = model
        self.settings = settings
        self.alpha = 1e-3
        self.expected_alpha = 1.0
        self.interpolation_variables_list = [KratosMultiphysics.KratosGlobals.GetVariable(i) for i in self.settings["interpolation_variables_list"].GetStringArray()]

    def ExecuteInitialize(self):
        self.origin_model_part_file_name = self.settings["origin_model_part_file_name"].GetString()
        destination_model_part_name = self.settings["destination_model_part_name"].GetString()
        self.destination_model_part = self.model.GetModelPart(destination_model_part_name)

        self.InterpolateVelocityWithMA()

    def InterpolateVelocityWithMA(self):
        #Import Origin Model Part 
        self.origin_model_part = self.model.CreateModelPart("OriginModelPart")
        for variable in self.destination_model_part.GetHistoricalVariablesNames():
            if KratosMultiphysics.KratosGlobals.HasVariable(variable):
                self.origin_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.KratosGlobals.GetVariable(variable))

        model_part_io = KratosMultiphysics.ModelPartIO(self.origin_model_part_file_name)
        model_part_io.ReadModelPart(self.origin_model_part) 

        #Set velocity field to origin model part
        for node in self.origin_model_part.Nodes:
            for variable in self.interpolation_variables_list:
                node.SetValue(variable, self.alpha * node.GetSolutionStepValue(variable))

        #Interpolate velocity to destination model part 
        interpolation = KratosMA.NodalValuesInterpolationProcess2D(self.origin_model_part, self.destination_model_part)
        interpolation.Execute()
        for node in self.destination_model_part.Nodes:
            for variable in self.interpolation_variables_list:
                node.SetSolutionStepValue(variable, node.GetValue(variable))

    def ExecuteInitializeSolutionStep(self):
        self.old_alpha = self.alpha
        self.alpha += 1e-3
        if self.alpha < self.expected_alpha:
            for node in self.destination_model_part.Nodes:
                for variable in self.interpolation_variables_list:
                    node.SetSolutionStepValue(variable, self.alpha * node.GetSolutionStepValue(variable) / self.old_alpha)
                    node.Fix(variable)