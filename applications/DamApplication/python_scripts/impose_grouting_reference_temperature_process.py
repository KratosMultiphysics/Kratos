import KratosMultiphysics
import KratosMultiphysics.DamApplication as KratosDam

def Factory(settings, Model):
    if(not isinstance(settings, KratosMultiphysics.Parameters)):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return ImposeGroutingReferenceTemperatureProcess(Model, settings["Parameters"])

class ImposeGroutingReferenceTemperatureProcess(KratosMultiphysics.Process):
    def __init__(self, Model, settings ):

        KratosMultiphysics.Process.__init__(self)
        model_part = Model[settings["model_part_name"].GetString()]

        self.process = KratosDam.DamGroutingReferenceTemperatureProcess(model_part, settings)


    def ExecuteBeforeSolutionLoop(self):
        self.process.ExecuteBeforeSolutionLoop()

    def ExecuteInitializeSolutionStep(self):
        pass

    def ExecuteFinalizeSolutionStep(self):
        self.process.ExecuteFinalizeSolutionStep()
