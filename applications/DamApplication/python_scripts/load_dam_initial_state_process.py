import KratosMultiphysics as KM
from KratosMultiphysics.DamApplication import dam_initial_state_file


def Factory(settings, model):
    if not isinstance(settings, KM.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return LoadDamInitialStateProcess(model, settings["Parameters"])


class LoadDamInitialStateProcess(KM.Process):
    """Loads selected fields from a lightweight Dam initial-state file."""

    def __init__(self, model, settings):
        super().__init__()
        settings.ValidateAndAssignDefaults(self._GetDefaultParameters())
        model_part_name = settings["model_part_name"].GetString()
        if model_part_name == "":
            raise Exception("LoadDamInitialStateProcess: 'model_part_name' is required.")
        if settings["input_file_name"].GetString() == "":
            raise Exception("LoadDamInitialStateProcess: 'input_file_name' is required.")
        self.model_part = model[model_part_name]
        self.input_file_name = settings["input_file_name"].GetString()
        self.transfer_initial_stresses = settings["transfer_initial_stresses"].GetBool()
        self.transfer_initial_temperature = settings["transfer_initial_temperature"].GetBool()
        self.transfer_reference_temperature = settings[
            "transfer_reference_temperature"].GetBool()
        self.reference_temperature_policy = settings[
            "reference_temperature_policy"].GetString()
        self.echo_level = settings["echo_level"].GetInt()
        if not (self.transfer_initial_stresses or self.transfer_initial_temperature or
                self.transfer_reference_temperature):
            raise Exception(
                "LoadDamInitialStateProcess: select at least one field to import.")
        if self.reference_temperature_policy not in ("preserve", "initial_temperature"):
            raise Exception(
                "LoadDamInitialStateProcess: unsupported reference_temperature_policy '{}'."
                .format(self.reference_temperature_policy))
        self._applied = False

        # Used by DamAnalysis to reject duplicate initial-state processors and an
        # accidental combination with the legacy previous-results mapping.
        self.transfers_initial_stress = self.transfer_initial_stresses
        self.transfers_initial_temperature = self.transfer_initial_temperature
        self.transfers_reference_temperature = (
            self.transfer_reference_temperature or
            self.reference_temperature_policy == "initial_temperature")

    @classmethod
    def _GetDefaultParameters(cls):
        return KM.Parameters("""{
            "model_part_name": "",
            "input_file_name": "",
            "transfer_initial_stresses": true,
            "transfer_initial_temperature": false,
            "transfer_reference_temperature": false,
            "reference_temperature_policy": "preserve",
            "echo_level": 0
        }""")

    def ExecuteBeforeSolutionLoop(self):
        if self._applied:
            raise RuntimeError(
                "LoadDamInitialStateProcess: the initial state was already applied. "
                "The process must run only once.")
        data = dam_initial_state_file.ReadFile(self.input_file_name)
        dam_initial_state_file.ValidateDestination(
            self.model_part,
            data,
            self.transfer_initial_stresses,
            self.transfer_initial_temperature,
            self.transfer_reference_temperature,
            self.reference_temperature_policy)
        dam_initial_state_file.ApplyTemperatures(
            self.model_part,
            data,
            self.transfer_initial_temperature,
            self.transfer_reference_temperature,
            self.reference_temperature_policy)
        if self.transfer_initial_stresses:
            dam_initial_state_file.ApplyStresses(self.model_part, data)
        self._applied = True
        if self.echo_level > 0:
            KM.Logger.PrintInfo(
                "LoadDamInitialStateProcess",
                "Loaded Dam initial-state file '{}'.".format(self.input_file_name))
            if self.transfer_initial_stresses:
                KM.Logger.PrintInfo(
                    "LoadDamInitialStateProcess",
                    "Imported stresses are interpreted as '{}'.".format(
                        data["stress_source"]))
