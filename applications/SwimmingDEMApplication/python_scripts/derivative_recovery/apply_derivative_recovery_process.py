import KratosMultiphysics
import KratosMultiphysics.SwimmingDEMApplication as SDEM

def Factory(settings, Model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return DerivativeRecoveryProcess(Model, settings["Parameters"])

class DerivativeRecoveryProcess(KratosMultiphysics.Process):

    def __init__(self, Model, settings):

        KratosMultiphysics.Process.__init__(self)

        default_parameters = KratosMultiphysics.Parameters( """
        {
            "model_part_name"                        : "FluidModelPart",
            "distance_factor"                        : 2.0,
            "distance_threshold"                     : 0.01,
            "continuous_distance"                    : true,
            "check_at_each_time_step"                : false,
            "avoid_almost_empty_elements"            : true,
            "deactivate_full_negative_elements"      : true,
            "recover_original_distance_at_each_step" : false
        }  """ )

        settings.ValidateAndAssignDefaults(default_parameters);

        self.fluid_model_part = Model[settings["model_part_name"].GetString()]
        self.DerivativeRecoveryProcess = SDEM.DerivativeRecoveryProcess(self.fluid_model_part, settings)


    def ExecuteInitialize(self):
        self.DerivativeRecoveryProcess.ExecuteInitialize()


    def ExecuteBeforeSolutionLoop(self):
        self.DerivativeRecoveryProcess.ExecuteBeforeSolutionLoop()


    def ExecuteInitializeSolutionStep(self):
        self.DerivativeRecoveryProcess.ExecuteInitializeSolutionStep()


    def ExecuteFinalizeSolutionStep(self):
        self.DerivativeRecoveryProcess.ExecuteFinalizeSolutionStep()
