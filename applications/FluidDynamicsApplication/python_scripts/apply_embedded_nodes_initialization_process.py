import KratosMultiphysics
import KratosMultiphysics.FluidDynamicsApplication as KratosFluid

def Factory(settings, Model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return ApplyEmbeddedNodesInitializationProcess(Model, settings["Parameters"])

class ApplyEmbeddedNodesInitializationProcess(KratosMultiphysics.Process):

    def __init__(self, Model, settings):

        KratosMultiphysics.Process.__init__(self)

        default_parameters = KratosMultiphysics.Parameters( """
        {
            "mesh_id"                   : 0,
            "model_part_name"           : "CHOOSE_FLUID_MODELPART_NAME",
            "max_iteration"             : 10
        }  """ )

        settings.ValidateAndAssignDefaults(default_parameters);

        self.fluid_model_part = Model[settings["model_part_name"].GetString()]
        self.max_iteration = settings["max_iteration"].GetInt()

        self.EmbeddedNodesInitializationProcess = KratosFluid.EmbeddedNodesInitializationProcess(self.fluid_model_part,
                                                                                                 self.max_iteration)


    def ExecuteInitializeSolutionStep(self):
        self.EmbeddedNodesInitializationProcess.ExecuteInitializeSolutionStep()
