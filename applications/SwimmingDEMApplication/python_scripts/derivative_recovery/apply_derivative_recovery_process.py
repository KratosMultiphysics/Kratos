import KratosMultiphysics
import KratosMultiphysics.SwimmingDEMApplication as SDEM

def Factory(settings, Model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return DerivativeRecoveryProcess(Model, settings["Parameters"])

class DerivativeRecoveryProcess(KratosMultiphysics.Process):
    @staticmethod
    def SDEMGetClassByName(name):
        if hasattr(SDEM, name):
            return getattr(SDEM, name)
        else:
            raise Exception('No python class or function called ' + name + 'found in SwimmingDEMApplication.')

    def __init__(self, Model, settings):

        KratosMultiphysics.Process.__init__(self)

        default_parameters = KratosMultiphysics.Parameters( """
        {
            "model_part_name" : "FluidModelPart",
            "recoverers" : [{}]
        }  """ )

        settings.ValidateAndAssignDefaults(default_parameters);

        self.fluid_model_part = Model[settings["model_part_name"].GetString()]
        self.list_of_recoverers = []

        for recoverer in settings["recoverers"].values():
            recovery_class = DerivativeRecoveryProcess.SDEMGetClassByName(recoverer["recoverer_name"].GetString())
            recoverer_instance = recovery_class(self.fluid_model_part, recoverer)
            self.list_of_recoverers.append(recoverer_instance)

    def ExecuteInitialize(self):
        for recoverer in self.list_of_recoverers:
            recoverer.ExecuteInitialize()

    def ExecuteBeforeSolutionLoop(self):
        for recoverer in self.list_of_recoverers:
            recoverer.ExecuteBeforeSolutionLoop()

    def ExecuteInitializeSolutionStep(self):
        for recoverer in self.list_of_recoverers:
            recoverer.ExecuteInitializeSolutionStep()

    def ExecuteFinalizeSolutionStep(self):
        for recoverer in self.list_of_recoverers:
            recoverer.ExecuteFinalizeSolutionStep()