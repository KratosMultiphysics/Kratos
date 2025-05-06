import KratosMultiphysics
import KratosMultiphysics.LaserDrillingApplication as KratosLaserDrilling


def Factory(settings, Model):
    if type(settings) != KratosMultiphysics.Parameters:  # TODO: Ruff E721
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return ElementDeactivationProcess(Model, settings["Parameters"])


class ElementDeactivationProcess(KratosMultiphysics.Process):
    """
    This process checks the thermal energy at each element and deactivates it if an input threshold is exceeded.
    """

    def __init__(self, Model, settings):
        KratosMultiphysics.Process.__init__(self)

        # Compare with default settings
        default_settings = KratosMultiphysics.Parameters(r"""{
            "model_part_name"          : "please_specify_model_part_name",
            "interval"                 : [0.0, 1e30],
            "thermal_energy_per_volume_threshold" : 1.0e20,
            "thermal_counter_threshold": 0,
            "alpha_threshold": 0.8,
            "decomposition_law" : "Prout-Tompkins",
            "A" : 8.72e6,
            "Ea": 7.96e4,
            "R" : 8.3145,
            "m" : 0.947,
            "n" : 3.692,
            "decomposition_law_reference_temperature": 400.0,
            "decomposition_law_constant_1": 1e-7,
            "decomposition_law_constant_2": 0.007
        }""")  # TODO: move to a "default parameters" file

        # Assign this here since it will change the "interval" prior to validation
        self.interval = KratosMultiphysics.IntervalUtility(settings)

        settings.ValidateAndAssignDefaults(default_settings)

        self.model_part = Model[settings["model_part_name"].GetString()]

        # Set parameters of the processes
        params = KratosMultiphysics.Parameters("{}")
        params.AddValue("model_part_name", settings["model_part_name"])
        params.AddValue("thermal_counter_threshold", settings["thermal_counter_threshold"])
        params.AddValue("thermal_energy_per_volume_threshold", settings["thermal_energy_per_volume_threshold"])
        params.AddValue("alpha_threshold", settings["alpha_threshold"])
        params.AddValue("decomposition_law", settings["decomposition_law"])
        params.AddValue("decomposition_law_reference_temperature", settings["decomposition_law_reference_temperature"])
        params.AddValue("decomposition_law_constant_1", settings["decomposition_law_constant_1"])
        params.AddValue("decomposition_law_constant_2", settings["decomposition_law_constant_2"])

        self.model_part.ProcessInfo[KratosLaserDrilling.THERMAL_ENERGY_PER_VOLUME_THRESHOLD] = 0.0
        self.model_part.ProcessInfo.SetValue(
            KratosLaserDrilling.THERMAL_ENERGY_PER_VOLUME_THRESHOLD,
            settings["thermal_energy_per_volume_threshold"].GetDouble(),
        )
        self.model_part.ProcessInfo[KratosLaserDrilling.ALPHA_THRESHOLD] = 0.0
        self.model_part.ProcessInfo.SetValue(
            KratosLaserDrilling.ALPHA_THRESHOLD, settings["alpha_threshold"].GetDouble()
        )
        self.model_part.ProcessInfo[KratosLaserDrilling.DECOMPOSITION_LAW] = ""
        self.model_part.ProcessInfo.SetValue(
            KratosLaserDrilling.DECOMPOSITION_LAW, settings["decomposition_law"].GetString()
        )
        # Prout-Tompkins law
        self.model_part.ProcessInfo[KratosLaserDrilling.A_COEFFICIENT] = 0.0
        self.model_part.ProcessInfo.SetValue(KratosLaserDrilling.A_COEFFICIENT, settings["A"].GetDouble())

        self.model_part.ProcessInfo[KratosLaserDrilling.EA_COEFFICIENT] = 0.0
        self.model_part.ProcessInfo.SetValue(KratosLaserDrilling.EA_COEFFICIENT, settings["Ea"].GetDouble())

        self.model_part.ProcessInfo[KratosLaserDrilling.R_GAS_CONSTANT] = 0.0
        self.model_part.ProcessInfo.SetValue(KratosLaserDrilling.R_GAS_CONSTANT, settings["R"].GetDouble())

        self.model_part.ProcessInfo[KratosLaserDrilling.M_COEFFICIENT] = 0.0
        self.model_part.ProcessInfo.SetValue(KratosLaserDrilling.M_COEFFICIENT, settings["m"].GetDouble())

        self.model_part.ProcessInfo[KratosLaserDrilling.N_COEFFICIENT] = 0.0
        self.model_part.ProcessInfo.SetValue(KratosLaserDrilling.N_COEFFICIENT, settings["n"].GetDouble())

        # Custom law
        self.model_part.ProcessInfo[KratosLaserDrilling.DECOMPOSITION_LAW_REFERENCE_TEMPERATURE] = 0.0
        self.model_part.ProcessInfo.SetValue(
            KratosLaserDrilling.DECOMPOSITION_LAW_REFERENCE_TEMPERATURE,
            settings["decomposition_law_reference_temperature"].GetDouble(),
        )

        self.model_part.ProcessInfo[KratosLaserDrilling.DECOMPOSITION_LAW_CONSTANT_1] = 0.0
        self.model_part.ProcessInfo.SetValue(
            KratosLaserDrilling.DECOMPOSITION_LAW_CONSTANT_1, settings["decomposition_law_constant_1"].GetDouble()
        )

        self.model_part.ProcessInfo[KratosLaserDrilling.DECOMPOSITION_LAW_CONSTANT_2] = 0.0
        self.model_part.ProcessInfo.SetValue(
            KratosLaserDrilling.DECOMPOSITION_LAW_CONSTANT_2, settings["decomposition_law_constant_2"].GetDouble()
        )

        # Set process
        self.process = KratosLaserDrilling.ElementDeactivationProcess(self.model_part, params)

    def ExecuteInitializeSolutionStep(self):
        current_time = self.model_part.ProcessInfo[KratosMultiphysics.TIME]

        if self.interval.IsInInterval(current_time):
            self.process.ExecuteInitializeSolutionStep()
