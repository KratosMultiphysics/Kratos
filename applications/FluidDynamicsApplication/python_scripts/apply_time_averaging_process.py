import KratosMultiphysics
import KratosMultiphysics.FluidDynamicsApplication as KratosFluid

def Factory(settings, Model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return ApplyTimeAveragingProcess(Model, settings["Parameters"])

class ApplyTimeAveragingProcess(KratosMultiphysics.Process):

    def __init__(self, Model, settings):

        KratosMultiphysics.Process.__init__(self)

        default_parameters = KratosMultiphysics.Parameters("""
        {
            "mesh_id"                   : 0,
            "model_part_name" : ""
        }
        """)

        settings.ValidateAndAssignDefaults(default_parameters);

        # Get model part name
        self.fluid_model_part = Model[settings["model_part_name"].GetString()]
        self.TimeAveragingProcess = KratosFluid.TimeAveragingProcess(self.fluid_model_part)



    def ExecuteInitialize(self):
        self.TimeAveragingProcess.ExecuteInitialize()

    def ExecuteFinalizeSolutionStep(self):
        self.TimeAveragingProcess.ExecuteFinalizeSolutionStep()


