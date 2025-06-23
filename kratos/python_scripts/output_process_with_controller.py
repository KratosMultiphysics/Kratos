# Importing the Kratos Library
import KratosMultiphysics

from KratosMultiphysics.model_parameters_factory import KratosModelParametersFactory


def Factory(settings, Model):
    if not isinstance(settings, KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return OutputProcessWithController(Model, settings["Parameters"])

class OutputProcessWithController(KratosMultiphysics.OutputProcess):
    """ This process divides its responsabilities between two objects:
        a controller that decides whether current time step must be
        printed or not and a print process that controls what will
        be printed and in which format.
    """
    def __init__(self, model, params):
        KratosMultiphysics.OutputProcess.__init__(self)
        default_settings = self.GetDefaultParameters()
        self.params = params
        self.model = model
        self.params.ValidateAndAssignDefaults(default_settings)
        factory = KratosModelParametersFactory(self.model)
        self.controller = factory.ConstructItem(self.params["controller_settings"])
        self.output_process = factory.ConstructItem(self.params["output_settings"])

    @classmethod
    def GetDefaultParameters(self):
        return KratosMultiphysics.Parameters('''{
            "controller_settings" : {},
            "output_settings" : {}
        }''')

    def Check(self):
        self.controller.Check()
        self.output_process.Check()

    def ExecuteInitialize(self):
        self.output_process.ExecuteInitialize()

    def ExecuteBeforeSolutionLoop(self):
        self.output_process.ExecuteBeforeSolutionLoop()

    def ExecuteInitializeSolutionStep(self):
        self.output_process.ExecuteInitializeSolutionStep()

    def ExecuteBeforeOutputStep(self):
        self.output_process.ExecuteBeforeOutputStep()

    def ExecuteAfterOutputStep(self):
        self.output_process.ExecuteAfterOutputStep()

    def ExecuteFinalizeSolutionStep(self):
        self.output_process.ExecuteFinalizeSolutionStep()

    def ExecuteFinalize(self):
        self.output_process.ExecuteFinalize()

    def IsOutputStep(self):
        self.controller.Evaluate()

    def PrintOutput(self):
        self.output_process.PrintOutput()
