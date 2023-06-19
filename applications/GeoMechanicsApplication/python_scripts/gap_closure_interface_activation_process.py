import KratosMultiphysics
import KratosMultiphysics.GeoMechanicsApplication as KratosGeo

def Factory(settings, Model):
    if not isinstance(settings, KratosMultiphysics.Parameters):
        raise TypeError("expected input shall be a Parameters object, encapsulating a json string")
    return GapClosureInterfaceActivationProcess(Model, settings["Parameters"])

## All the python processes should be derived from "python_process"

class GapClosureInterfaceActivationProcess(KratosMultiphysics.Process):
    def __init__(self, Model, settings ):
        KratosMultiphysics.Process.__init__(self)

        model_part = Model[settings["model_part_name"].GetString()]

        params = KratosMultiphysics.Parameters("{}")
        params.AddValue("model_part_name",settings["model_part_name"])
        params.AddValue("gap_width_threshold",settings["gap_width_threshold"])
        params.AddValue("consider_gap_closure",settings["consider_gap_closure"])
        self.process = KratosGeo.GapClosureInterfaceProcess(model_part, params)

    def ExecuteInitialize(self):
        self.process.ExecuteInitialize()

    def ExecuteInitializeSolutionStep(self):
        self.process.ExecuteInitializeSolutionStep()

    def ExecuteFinalizeSolutionStep(self):
        self.process.ExecuteFinalizeSolutionStep()

    def ExecuteFinalize(self):
        self.process.ExecuteFinalize()
