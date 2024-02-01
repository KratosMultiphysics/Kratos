import KratosMultiphysics
import KratosMultiphysics.FluidDynamicsApplication as KratosFluid

def Factory(settings, Model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    if settings["Parameters"].Has("mesh_id"):
        settings["Parameters"].RemoveValue("mesh_id")
        KratosMultiphysics.Logger.PrintWarning("ApplyEmbeddedPostprocessProcess", "mesh_id is a legacy setting. Please remove mesh_id from your parameters")
    return ApplyEmbeddedPostprocessProcess(Model, settings["Parameters"])

class ApplyEmbeddedPostprocessProcess(KratosMultiphysics.Process):

    def __init__(self, Model, settings):

        KratosMultiphysics.Process.__init__(self)

        default_parameters = KratosMultiphysics.Parameters( """
        {
            "model_part_name"           : "CHOOSE_FLUID_MODELPART_NAME"
        }  """ )

        settings.ValidateAndAssignDefaults(default_parameters);

        self.fluid_model_part = Model[settings["model_part_name"].GetString()]

        self.EmbeddedPostprocessProcess = KratosFluid.EmbeddedPostprocessProcess(self.fluid_model_part)


    def ExecuteFinalizeSolutionStep(self):
        self.EmbeddedPostprocessProcess.ExecuteFinalizeSolutionStep()
