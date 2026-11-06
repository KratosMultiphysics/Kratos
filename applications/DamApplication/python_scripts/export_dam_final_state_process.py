import KratosMultiphysics as KM
from KratosMultiphysics.DamApplication import dam_initial_state_file


def Factory(settings, model):
    if not isinstance(settings, KM.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return ExportDamFinalStateProcess(model, settings["Parameters"])


class ExportDamFinalStateProcess(KM.Process):
    """Writes a Dam final-state file at the end of the analysis."""

    def __init__(self, model, settings):
        super().__init__()
        settings.ValidateAndAssignDefaults(self._GetDefaultParameters())
        model_part_name = settings["model_part_name"].GetString()
        if model_part_name == "":
            raise Exception(
                "ExportDamFinalStateProcess: 'model_part_name' is required.")
        if settings["output_file_name"].GetString() == "":
            raise Exception(
                "ExportDamFinalStateProcess: 'output_file_name' is required.")
        self.model_part = model[model_part_name]
        self.export_stresses = settings["export_stresses"].GetBool()
        self.export_temperature = settings["export_temperature"].GetBool()
        self.export_reference_temperature = settings[
            "export_reference_temperature"].GetBool()
        self.output_file_name = settings["output_file_name"].GetString()
        self.echo_level = settings["echo_level"].GetInt()
        if not (self.export_stresses or self.export_temperature or
                self.export_reference_temperature):
            raise Exception(
                "ExportDamFinalStateProcess: select at least one field to export.")
        self._exported = False

    @classmethod
    def _GetDefaultParameters(cls):
        return KM.Parameters("""{
            "model_part_name": "",
            "export_stresses": true,
            "export_temperature": false,
            "export_reference_temperature": false,
            "output_file_name": "",
            "echo_level": 0
        }""")

    def ExecuteFinalize(self):
        if self._exported:
            return
        data = dam_initial_state_file.BuildFileData(
            self.model_part,
            self.export_stresses,
            self.export_temperature,
            self.export_reference_temperature)
        dam_initial_state_file.WriteFile(self.output_file_name, data)
        self._exported = True
        if self.echo_level > 0:
            KM.Logger.PrintInfo(
                "ExportDamFinalStateProcess",
                "Wrote Dam final-state file '{}'.".format(self.output_file_name))
