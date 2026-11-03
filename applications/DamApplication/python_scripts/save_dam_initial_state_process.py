import KratosMultiphysics as KM
from KratosMultiphysics.DamApplication import dam_initial_state_file


def Factory(settings, model):
    if not isinstance(settings, KM.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return SaveDamInitialStateProcess(model, settings["Parameters"])


class SaveDamInitialStateProcess(KM.Process):
    """Writes a lightweight Dam initial-state file at the end of the analysis."""

    def __init__(self, model, settings):
        super().__init__()
        settings.ValidateAndAssignDefaults(self._GetDefaultParameters())
        model_part_name = settings["model_part_name"].GetString()
        if model_part_name == "":
            raise Exception("SaveDamInitialStateProcess: 'model_part_name' is required.")
        if settings["output_file_name"].GetString() == "":
            raise Exception("SaveDamInitialStateProcess: 'output_file_name' is required.")
        self.model_part = model[model_part_name]
        self.save_initial_stresses = settings["save_initial_stresses"].GetBool()
        self.save_initial_temperature = settings["save_initial_temperature"].GetBool()
        self.save_reference_temperature = settings["save_reference_temperature"].GetBool()
        self.stress_source = settings["stress_source"].GetString()
        self.output_file_name = settings["output_file_name"].GetString()
        self.echo_level = settings["echo_level"].GetInt()
        if not (self.save_initial_stresses or self.save_initial_temperature or
                self.save_reference_temperature):
            raise Exception(
                "SaveDamInitialStateProcess: select at least one field to save.")
        if self.save_initial_stresses and self.stress_source not in (
                "current_cauchy", "stored_initial_stress"):
            raise Exception(
                "SaveDamInitialStateProcess: unsupported stress_source '{}'.".format(
                    self.stress_source))
        self._saved = False

    @classmethod
    def _GetDefaultParameters(cls):
        return KM.Parameters("""{
            "model_part_name": "",
            "save_initial_stresses": true,
            "save_initial_temperature": false,
            "save_reference_temperature": false,
            "stress_source": "current_cauchy",
            "output_file_name": "",
            "echo_level": 0
        }""")

    def ExecuteFinalize(self):
        if self._saved:
            return
        data = dam_initial_state_file.BuildFileData(
            self.model_part,
            self.save_initial_stresses,
            self.save_initial_temperature,
            self.save_reference_temperature,
            self.stress_source)
        dam_initial_state_file.WriteFile(self.output_file_name, data)
        self._saved = True
        if self.echo_level > 0:
            KM.Logger.PrintInfo(
                "SaveDamInitialStateProcess",
                "Wrote Dam initial-state file '{}'.".format(self.output_file_name))
