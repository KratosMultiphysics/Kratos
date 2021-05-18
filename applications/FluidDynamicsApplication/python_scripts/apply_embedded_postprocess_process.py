import KratosMultiphysics
import KratosMultiphysics.FluidDynamicsApplication as KratosFluid

def Factory(settings, Model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return ApplyEmbeddedPostprocessProcess(Model, settings["Parameters"])

class ApplyEmbeddedPostprocessProcess(KratosMultiphysics.Process):

    def __init__(self, Model, settings):

        KratosMultiphysics.Process.__init__(self)

        default_parameters = KratosMultiphysics.Parameters( """
        {
            "mesh_id"                   : 0,
            "model_part_name"           : "CHOOSE_FLUID_MODELPART_NAME"
        }  """ )

        settings.ValidateAndAssignDefaults(default_parameters);

        self.fluid_model_part = Model[settings["model_part_name"].GetString()]

        self.EmbeddedPostprocessProcess = KratosFluid.EmbeddedPostprocessProcess(self.fluid_model_part)


    def ExecuteFinalizeSolutionStep(self):
        self.EmbeddedPostprocessProcess.ExecuteFinalizeSolutionStep()
