
import KratosMultiphysics
import KratosMultiphysics.GeoMechanicsApplication as KratosGeo

def Factory(settings, model):
    if not isinstance(settings, KratosMultiphysics.Parameters):
        raise TypeError("expected input shall be a Parameters object, encapsulating a json string")
    return ApplyCPhiReductionProcess(model, settings["Parameters"])

## All the python processes should be derived from "python_process"

class ApplyCPhiReductionProcess(KratosMultiphysics.Process):
    def __init__(self, model, settings ):
        KratosMultiphysics.Process.__init__(self)

        model_part = model[settings["model_part_name"].GetString()]

        params = KratosMultiphysics.Parameters("{}")
        params.AddValue("model_part_name",settings["model_part_name"])

        self.process = KratosGeo.ApplyCPhiReductionProcess(model_part,params)

    def ExecuteInitializeSolutionStep(self):
        self.process.ExecuteInitializeSolutionStep()

    def ExecuteFinalize(self):
        self.process.ExecuteFinalize()
