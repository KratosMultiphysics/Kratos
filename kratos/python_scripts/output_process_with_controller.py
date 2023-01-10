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
        self.print_process = factory.ConstructItem(self.params["output_settings"])

    @classmethod
    def GetDefaultParameters(self):
        return KratosMultiphysics.Parameters('''{
            "controller_settings" : {},
            "output_settings" : {}
        }''')

    def ExecuteInitialize(self):
        self.print_process.ExecuteInitialize()

    def ExecuteBeforeSolutionLoop(self):
        self.print_process.ExecuteBeforeSolutionLoop()

    def ExecuteInitializeSolutionStep(self):
        self.print_process.ExecuteInitializeSolutionStep()

    def ExecuteBeforeOutputStep(self):
        self.print_process.ExecuteBeforeOutputStep()

    def ExecuteAfterOutputStep(self):
        self.print_process.ExecuteAfterOutputStep()

    def ExecuteFinalizeSolutionStep(self):
        self.print_process.ExecuteFinalizeSolutionStep()

    def ExecuteFinalize(self):
        self.print_process.ExecuteFinalize()

    def IsOutputStep(self):
        self.controller.IsOutputStep()

    def PrintOutput(self):
        self.print_process.PrintOutput()
