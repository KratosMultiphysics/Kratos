import KratosMultiphysics
import KratosMultiphysics.FluidDynamicsApplication as KratosFluid

def Factory(settings, Model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return ApplyEmbeddedSkinVisualizationProcess(Model, settings["Parameters"])

class ApplyEmbeddedSkinVisualizationProcess(KratosMultiphysics.Process):

    def __init__(self, Model, settings):

        KratosMultiphysics.Process.__init__(self)

        default_parameters = KratosMultiphysics.Parameters( """
        {
            "model_part_name"                     : "origin_model_part",
            "visualization_model_part_name"       : "visualization_model_part",
            "shape_functions"                     : "standard",
            "reform_model_part_at_each_time_step" : false,
            "visualization_variables"             : ["VELOCITY","PRESSURE"]
        }  """ )

        settings.ValidateAndAssignDefaults(default_parameters);

        self.origin_model_part = Model[settings["model_part_name"].GetString()]
        self.visualization_model_part = ModelPart(settings["visualization_model_part_name"].GetString())

        aux_params = settings
        aux_params.RemoveValue("model_part_name")
        aux_params.RemoveValue("visualization_model_part_name")

        self.EmbeddedSkinVisualizationProcess = KratosFluid.EmbeddedSkinVisualizationProcess(
            self.origin_model_part,
            self.visualization_model_part,
            aux_params)

    def ExecuteInitialize(self):
        self.EmbeddedSkinVisualizationProcess.ExecuteInitialize()

    def ExecuteInitializeSolutionStep(self):
        self.EmbeddedSkinVisualizationProcess.ExecuteInitializeSolutionStep()

    def ExecuteBeforeOutputStep(self):
        self.EmbeddedSkinVisualizationProcess.ExecuteBeforeOutputStep()

    def ExecuteFinalizeSolutionStep(self):
        self.EmbeddedSkinVisualizationProcess.ExecuteFinalizeSolutionStep()
