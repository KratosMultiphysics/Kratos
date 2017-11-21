import KratosMultiphysics
import KratosMultiphysics.FluidDynamicsApplication as KratosFluid
import math

def Factory(settings, Model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return ApplyDistanceModificationProcess(Model, settings["Parameters"])

class ApplyDistanceModificationProcess(KratosMultiphysics.Process):

    def __init__(self, Model, settings):

        KratosMultiphysics.Process.__init__(self)

        default_parameters = KratosMultiphysics.Parameters( """
        {
            "mesh_id"                                   : 0,
            "model_part_name"                           : "CHOOSE_FLUID_MODELPART_NAME",
            "check_at_each_time_step"                   : false,
            "deactivate_full_negative_elements"         : true,
            "recover_original_distance_at_each_step"    : false
        }  """ )

        settings.ValidateAndAssignDefaults(default_parameters);

        self.fluid_model_part = Model[settings["model_part_name"].GetString()]
        self.DistanceModificationProcess = KratosFluid.DistanceModificationProcess(self.fluid_model_part, settings)


    def ExecuteInitialize(self):
        self.DistanceModificationProcess.ExecuteInitialize()


    def ExecuteBeforeSolutionLoop(self):
        self.DistanceModificationProcess.ExecuteBeforeSolutionLoop()


    def ExecuteInitializeSolutionStep(self):
        self.DistanceModificationProcess.ExecuteInitializeSolutionStep()


    def ExecuteFinalizeSolutionStep(self):
        self.DistanceModificationProcess.ExecuteFinalizeSolutionStep()
