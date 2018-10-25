import KratosMultiphysics
import KratosMultiphysics.FluidDynamicsApplication as KratosFluid
import math

def Factory(settings, Model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return ApplyMassConservationCheckProcess(Model, settings["Parameters"])

class ApplyMassConservationCheckProcess(KratosMultiphysics.Process):

    def __init__(self, Model, settings):

        KratosMultiphysics.Process.__init__(self)

        default_parameters = KratosMultiphysics.Parameters( """
        {
            "model_part_name"                        : "default_model_part_name",
            "is_switched_on"                         : true,
            "mass_computation_frequency"             : 20,
            "compare_to_initial_values"              : true,
            "write_to_log_file"                      : true
        }  """ )

        settings.ValidateAndAssignDefaults(default_parameters);

        self.fluid_model_part = Model[settings["model_part_name"].GetString()]
        self.MassConservationCheckProcess = KratosFluid.MassConservationCheckProcess(self.fluid_model_part, settings)


    def ExecuteInitialize(self):
        self.MassConservationCheckProcess.ExecuteInitialize()


    def ExecuteBeforeSolutionLoop(self):
        self.MassConservationCheckProcess.ExecuteBeforeSolutionLoop()


    def ExecuteInitializeSolutionStep(self):
        self.MassConservationCheckProcess.ExecuteInitializeSolutionStep()


    def ExecuteFinalizeSolutionStep(self):
        self.MassConservationCheckProcess.ExecuteFinalizeSolutionStep()
