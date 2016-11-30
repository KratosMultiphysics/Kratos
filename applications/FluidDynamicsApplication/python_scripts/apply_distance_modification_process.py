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
            "mesh_id"                   : 0,
            "model_part_name"           : "CHOOSE_FLUID_MODELPART_NAME",
            "check_at_each_time_step"   : false
        }  """ )

        settings.ValidateAndAssignDefaults(default_parameters);
        
        self.fluid_model_part = Model[settings["model_part_name"].GetString()]
        self.check_at_each_time_step = settings["check_at_each_time_step"].GetBool()
        
        self.DistanceModificationProcess = KratosFluid.DistanceModificationProcess(self.fluid_model_part,
                                                                                   self.check_at_each_time_step)


    def ExecuteInitialize(self):
        self.DistanceModificationProcess.ExecuteInitialize()
        
    
    def ExecuteBeforeSolutionLoop(self):
        self.DistanceModificationProcess.ExecuteBeforeSolutionLoop()
        
        
    def ExecuteInitializeSolutionStep(self):
        self.DistanceModificationProcess.ExecuteInitializeSolutionStep()
        
        
    def ExecuteFinalize(self):
        self.DistanceModificationProcess.ExecuteFinalize()
