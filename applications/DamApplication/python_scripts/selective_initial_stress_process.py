import KratosMultiphysics as KM
import KratosMultiphysics.DamApplication as KratosDam
from KratosMultiphysics.restart_utility import RestartUtility


def Factory(settings, model):
    if not isinstance(settings, KM.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return SelectiveInitialStressProcess(model, settings["Parameters"])


class SelectiveInitialStressProcess(KM.Process):
    """Applies selected Gauss-point stresses from a restart checkpoint to a new analysis.

    The source is loaded into an isolated Model and is never initialized, solved or
    modified. Only the source Cauchy stresses cross the model boundary. The
    destination keeps its own geometry, displacement datum, time, step, loads,
    boundary conditions and temperature fields.
    """

    def __init__(self, model, settings):
        super().__init__()
        settings.ValidateAndAssignDefaults(self._GetDefaultParameters())

        model_part_name = settings["model_part_name"].GetString()
        if model_part_name == "":
            raise Exception("SelectiveInitialStressProcess: 'model_part_name' is required.")

        restart_settings = settings["restart_settings"]
        if (not restart_settings.Has("input_filename") or
                restart_settings["input_filename"].GetString() == ""):
            raise Exception(
                "SelectiveInitialStressProcess: 'restart_settings.input_filename' is required.")

        self.destination_model_part = model[model_part_name]
        self.source_model_part_name = settings["source_model_part_name"].GetString()
        self.source_root_model_part_name = settings["source_root_model_part_name"].GetString()
        if self.source_root_model_part_name == "":
            raise Exception(
                "SelectiveInitialStressProcess: 'source_root_model_part_name' is required.")
        self.restart_settings = restart_settings
        self._applied = False

        # Used by DamAnalysis to reject an accidental combination with the legacy
        # previous-results displacement/stress mapping.
        self.transfers_initial_stress = True

    @classmethod
    def _GetDefaultParameters(cls):
        return KM.Parameters("""{
            "model_part_name": "",
            "source_model_part_name": "",
            "source_root_model_part_name": "MainModelPart",
            "restart_settings": {}
        }""")

    def ExecuteBeforeSolutionLoop(self):
        if self._applied:
            raise RuntimeError(
                "SelectiveInitialStressProcess: the initial stress state was already applied. "
                "The process must run only once.")
        self._applied = True

        source_model = KM.Model()
        try:
            source_root = source_model.CreateModelPart(self.source_root_model_part_name)
            RestartUtility(source_root, self.restart_settings).LoadRestart()

            if self.source_model_part_name == "":
                source_model_part = source_root
            else:
                if not source_root.HasSubModelPart(self.source_model_part_name):
                    raise RuntimeError(
                        "SelectiveInitialStressProcess: source model part '{}' was not found "
                        "in the restart checkpoint.".format(self.source_model_part_name))
                source_model_part = source_root.GetSubModelPart(self.source_model_part_name)

            KratosDam.SelectiveInitialStressUtility.AssignInitialStresses(
                source_model_part, self.destination_model_part)
        finally:
            # The C++ utility copies the selected stresses into the destination
            # InitialState objects. The source can be released unconditionally.
            del source_model
