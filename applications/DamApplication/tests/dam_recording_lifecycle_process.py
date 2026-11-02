import KratosMultiphysics


LAST_PROCESS = None


def Factory(settings, Model):
    global LAST_PROCESS
    LAST_PROCESS = RecordingLifecycleProcess(Model, settings["Parameters"])
    return LAST_PROCESS


class RecordingLifecycleProcess(KratosMultiphysics.Process):
    """Test-only process that records the lifecycle callbacks it receives."""

    def __init__(self, model, settings):
        super().__init__()
        self.settings = settings
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

    def ExecuteAfterOutputStep(self):
        self._record("ExecuteAfterOutputStep")

    def ExecuteFinalize(self):
        self._record("ExecuteFinalize")
