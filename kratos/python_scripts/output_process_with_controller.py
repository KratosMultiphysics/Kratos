# Importing the Kratos Library
import KratosMultiphysics


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
        controller_name = self.params["controller_settings"]["name"].GetString()
        controller_prototype = KratosMultiphysics.Registry[f"{controller_name}.Prototype"]
        output_process_name = self.params["output_settings"]["name"].GetString()
        output_settings_prototype = KratosMultiphysics.Registry[f"{output_process_name}.Prototype"]

        self.controller = controller_prototype.Create(self.model, self.params["controller_settings"]["parameters"])
        self.print_process = output_settings_prototype.Create(self.model, self.params["output_settings"]["parameters"])

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
