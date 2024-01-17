import KratosMultiphysics
import KratosMultiphysics.ConvectionDiffusionApplication as KratosConvDiff

def Factory(settings, Model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return ElementDeactivationProcess(Model, settings["Parameters"])


class ElementDeactivationProcess(KratosMultiphysics.Process):
    '''
    This process checks the thermal energy at each element and deactivates it if an input threshold is exceeded.
    '''

    def __init__(self, Model, settings):

        KratosMultiphysics.Process.__init__(self)

        # Compare with default settings
        default_settings = KratosMultiphysics.Parameters(r'''{
            "model_part_name"          : "please_specify_model_part_name",
            "interval"                 : [0.0, 1e30],
            "thermal_energy_per_volume_threshold" : 1.0e20,
            "thermal_counter_threshold": 0
        }''')

        # Assign this here since it will change the "interval" prior to validation
        self.interval = KratosMultiphysics.IntervalUtility(settings)

        settings.ValidateAndAssignDefaults(default_settings)

        self.model_part = Model[settings["model_part_name"].GetString()]

        # Set parameters of the processes
        params = KratosMultiphysics.Parameters("{}")
        params.AddValue("model_part_name", settings["model_part_name"])
        params.AddValue("thermal_energy_per_volume_threshold", settings["thermal_energy_per_volume_threshold"])
        params.AddValue("thermal_counter_threshold", settings["thermal_counter_threshold"])

        # Set process
        self.process = KratosConvDiff.ElementDeactivationProcess(self.model_part, params)

    def ExecuteInitializeSolutionStep(self):

        current_time = self.model_part.ProcessInfo[KratosMultiphysics.TIME]

        if self.interval.IsInInterval(current_time):
            self.process.ExecuteInitializeSolutionStep()

