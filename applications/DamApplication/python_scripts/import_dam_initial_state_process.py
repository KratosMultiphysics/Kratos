import KratosMultiphysics as KM
from KratosMultiphysics.DamApplication import dam_initial_state_file


def Factory(settings, model):
    if not isinstance(settings, KM.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return ImportDamInitialStateProcess(model, settings["Parameters"])


class ImportDamInitialStateProcess(KM.Process):
    """Loads selected fields from a Dam final-state file as an initial state."""

    def __init__(self, model, settings):
        super().__init__()
        settings.ValidateAndAssignDefaults(self._GetDefaultParameters())
        model_part_name = settings["model_part_name"].GetString()
        if model_part_name == "":
            raise Exception(
                "ImportDamInitialStateProcess: 'model_part_name' is required.")
        if settings["input_file_name"].GetString() == "":
            raise Exception(
                "ImportDamInitialStateProcess: 'input_file_name' is required.")
        self.model_part = model[model_part_name]
        self.input_file_name = settings["input_file_name"].GetString()
        self.import_stresses = settings["import_stresses"].GetBool()
        self.import_temperature = settings["import_temperature"].GetBool()
        self.import_reference_temperature = settings[
            "import_reference_temperature"].GetBool()
        self.reference_temperature_policy = settings[
            "reference_temperature_policy"].GetString()
        self.transfer_method = settings["transfer_method"].GetString()
        self.echo_level = settings["echo_level"].GetInt()
        if not (self.import_stresses or self.import_temperature or
                self.import_reference_temperature):
            raise Exception(
                "ImportDamInitialStateProcess: select at least one field to import.")
        if self.reference_temperature_policy not in ("preserve", "initial_temperature"):
            raise Exception(
                "ImportDamInitialStateProcess: unsupported reference_temperature_policy "
                "'{}'.".format(self.reference_temperature_policy))
        if self.transfer_method not in ("matching_mesh", "spatial_mapping"):
            raise Exception(
                "ImportDamInitialStateProcess: unsupported transfer_method "
                "'{}'.".format(self.transfer_method))
        self._applied = False

        # Used by the analysis to detect reference-field imports and to reject
        # combination with the legacy previous-results mapping.
        self.imports_reference_temperature = (
            self.import_reference_temperature or
            self.reference_temperature_policy == "initial_temperature")

    @classmethod
    def _GetDefaultParameters(cls):
        return KM.Parameters("""{
            "model_part_name": "",
            "input_file_name": "",
            "import_stresses": true,
            "import_temperature": false,
            "import_reference_temperature": false,
            "reference_temperature_policy": "preserve",
            "transfer_method": "matching_mesh",
            "echo_level": 0
        }""")

    def PrevalidateInitialState(self):
        """Reads and validates the file without modifying the destination.

        DamAnalysis calls this for every configured import before any process
        applies fields, so overlapping element/node assignments can be rejected
        globally before the model is modified.
        """
        if not hasattr(self, "_data"):
            data = dam_initial_state_file.ReadFile(self.input_file_name)
            if self.transfer_method == "spatial_mapping":
                from KratosMultiphysics.DamApplication import dam_initial_state_mapping
                dam_initial_state_mapping.ValidateSpatialDestination(
                    self.model_part,
                    data,
                    self.import_stresses,
                    self.import_temperature,
                    self.import_reference_temperature,
                    self.reference_temperature_policy)
                data = dam_initial_state_mapping.BuildSpatialMappedData(
                    self.model_part,
                    data,
                    self.import_stresses,
                    self.import_temperature,
                    self.import_reference_temperature,
                    self.reference_temperature_policy)
            else:
                dam_initial_state_file.ValidateDestination(
                    self.model_part,
                    data,
                    self.import_stresses,
                    self.import_temperature,
                    self.import_reference_temperature,
                    self.reference_temperature_policy)
            self._data = data
        return self._data

    def ExecuteBeforeSolutionLoop(self):
        if self._applied:
            raise RuntimeError(
                "ImportDamInitialStateProcess: the initial state was already applied. "
                "The process must run only once.")
        data = self.PrevalidateInitialState()
        dam_initial_state_file.ApplyTemperatures(
            self.model_part,
            data,
            self.import_temperature,
            self.import_reference_temperature,
            self.reference_temperature_policy)
        if self.import_stresses:
            dam_initial_state_file.ApplyStresses(self.model_part, data)
        self._applied = True
        if self.echo_level > 0:
            KM.Logger.PrintInfo(
                "ImportDamInitialStateProcess",
                "Loaded Dam final-state file '{}'.".format(self.input_file_name))
