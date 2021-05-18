import KratosMultiphysics
import KratosMultiphysics.DEMApplication as Dem

def Factory(settings, Model):
    if not isinstance(settings, KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return AutomaticDTProcess(Model, settings["Parameters"])

## All the processes python should be derived from "Process"

class AutomaticDTProcess(KratosMultiphysics.Process):
    def __init__(self, Model, settings ):
        KratosMultiphysics.Process.__init__(self)

        self.model_part = Model[settings["model_part_name"].GetString()]
        self.automatic_dt_process = Dem.AutomaticDTProcess(self.model_part, settings)

    def ExecuteBeforeSolutionLoop(self):
        self.automatic_dt_process.ExecuteBeforeSolutionLoop()