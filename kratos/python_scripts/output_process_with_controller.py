# Importing the Kratos Library
import KratosMultiphysics

# other imports
from importlib import import_module

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
        self.controller = self._CreateInstance(self.params["controller_settings"])
        self.print_process = self._CreateInstance(self.params["print_process_settings"])

    @classmethod
    def GetDefaultParameters(self):
        return KratosMultiphysics.Parameters('''{
            "controller_settings" : {},
            "print_process_settings" : {}
        }''')

    def _CreateInstance(self, settings):
        if not settings.Has("python_module"):
            raise NameError('"python_module" must be defined in your parameters. Check ams_output_process settings')
        python_module_name = settings["python_module"].GetString()

        if not settings.Has("kratos_module"):
            raise NameError('"kratos_module" must be defined in your parameters. Check ams_output_process settings')

        if not settings.Has("Parameters"):
            settings.AddEmptyValue("Parameters")

        kratos_module_name = settings["kratos_module"].GetString()
        if not kratos_module_name.startswith("KratosMultiphysics"):
            kratos_module_name = "KratosMultiphysics." + kratos_module_name

        full_module_name = kratos_module_name + "." + python_module_name
        python_module = import_module(full_module_name)

        return python_module.Factory(settings, self.model)

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
