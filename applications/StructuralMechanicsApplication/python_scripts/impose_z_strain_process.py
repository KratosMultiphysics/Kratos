import KratosMultiphysics
import KratosMultiphysics.StructuralMechanicsApplication as SMA

def Factory(settings, Model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return ImposeZStrainProcess(Model, settings["Parameters"])

## All the processes python should be derived from "Process"

class ImposeZStrainProcess(KratosMultiphysics.Process):
    def __init__(self, Model, settings ):
        KratosMultiphysics.Process.__init__(self)

        self.model_part = Model[settings["model_part_name"].GetString()]

        self.params = KratosMultiphysics.Parameters("{}")
        self.params.AddValue("model_part_name",settings["model_part_name"])
        self.params.AddValue("z_strain_value",settings["z_strain_value"])

        self.process = SMA.ImposeZStrainProcess(self.model_part, self.params)

    def ExecuteInitialize(self):

        self.process.ExecuteInitialize()

    def ExecuteInitializeSolutionStep(self):

        self.process.ExecuteInitializeSolutionStep()