"""Investigation-only adapter, executed once after destination solver initialization."""
import json

import KratosMultiphysics as KM
from KratosMultiphysics.restart_utility import RestartUtility
from test_dam_selective_initialization import snapshot


def Factory(settings, model):
    return SelectiveStressProbe(model, settings["Parameters"])


class SelectiveStressProbe(KM.Process):
    def __init__(self, model, settings):
        super().__init__()
        self.target = model[settings["model_part_name"].GetString()]
        self.settings = settings

    def ExecuteBeforeSolutionLoop(self):
        from dam_selective_state_probe import assign_stresses
        source_model = KM.Model()
        source = source_model.CreateModelPart("MainModelPart")
        RestartUtility(source, self.settings["restart"].Clone()).LoadRestart()
        source_part = source.GetSubModelPart("Parts_Concrete")
        before = snapshot(source_part)
        assign_stresses(source_part, self.target)
        after = snapshot(source_part)
        if before != after:
            raise RuntimeError("Selective extraction changed source state")
        with open("loaded_source.json", "w") as result:
            json.dump(before, result)
