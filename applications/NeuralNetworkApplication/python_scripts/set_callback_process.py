import KratosMultiphysics as KM
import KratosMultiphysics.NeuralNetworkApplication as NeuralNetwork
from KratosMultiphysics.NeuralNetworkApplication.neural_network_process import NeuralNetworkProcess


from importlib import import_module
from tensorflow.keras import layers

def Factory(settings):
    if not isinstance(settings, KM.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return SetCallbackProcess(settings["parameters"])

class SetCallbackProcess(NeuralNetworkProcess):

    def __init__(self, settings):
        super().__init__()
        """ The default constructor of the class

        Keyword arguments:
        self -- It signifies an instance of a class.
        model -- the container of the different model parts.
        settings -- Kratos parameters containing process settings.
        """

        default_settings = KM.Parameters("""{
            "help"                        : "This process adds a callback to a process.",
            "callback_type"               : "",
            "callback_parameters"         : {}
        }""")

        settings.ValidateAndAssignDefaults(default_settings)

        self.callback_type = settings["callback_type"].GetString()
        self.callback_parameters = settings["callback_parameters"]
        # Import layer classes
        module_name = "KratosMultiphysics.NeuralNetworkApplication." + self.callback_type
        class_module = import_module(module_name)
        callback_class = class_module.Factory(self.callback_parameters)

        # Build layer
        self.callback = callback_class.Build()

    def Callback(self):
        """ Processes to set callbacks to a network. """
        return self.callback