import KratosMultiphysics as KM


def Factory(settings, model):
    return SetInitialStateBeforeOutput(model, settings["Parameters"])


class SetInitialStateBeforeOutput(KM.Process):
    """Test adapter: apply the core process once after solver initialization."""

    def __init__(self, model, settings):
        super().__init__()
        self.model_part = model[settings["model_part_name"].GetString()]
        stress = KM.Vector(settings["initial_stress"].GetVector())
        self.process = KM.SetInitialStateProcess2D(self.model_part, stress, 1)

    def ExecuteBeforeSolutionLoop(self):
        self.process.ExecuteInitializeSolutionStep()
