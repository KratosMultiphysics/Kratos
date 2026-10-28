import KratosMultiphysics


class RecordingOutputProcess(KratosMultiphysics.OutputProcess):
    """Test-only output process that records the lifecycle calls it receives."""

    def __init__(self, model, settings):
        super().__init__()
        self.settings = settings
        self.settings.ValidateAndAssignDefaults(KratosMultiphysics.Parameters('{"record_name": ""}'))
        self.calls = []

    def _record(self, call_name):
        self.calls.append(call_name)

    def ExecuteInitialize(self):
        self._record("ExecuteInitialize")

    def ExecuteBeforeSolutionLoop(self):
        self._record("ExecuteBeforeSolutionLoop")

    def ExecuteInitializeSolutionStep(self):
        self._record("ExecuteInitializeSolutionStep")

    def ExecuteFinalizeSolutionStep(self):
        self._record("ExecuteFinalizeSolutionStep")

    def ExecuteBeforeOutputStep(self):
        self._record("ExecuteBeforeOutputStep")

    def IsOutputStep(self):
        self._record("IsOutputStep")
        return True

    def PrintOutput(self):
        self._record("PrintOutput")

    def ExecuteAfterOutputStep(self):
        self._record("ExecuteAfterOutputStep")

    def ExecuteFinalize(self):
        self._record("ExecuteFinalize")


def Factory(settings, model):
    return RecordingOutputProcess(model, settings["Parameters"])
