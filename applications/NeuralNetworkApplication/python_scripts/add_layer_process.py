import KratosMultiphysics as KM
import KratosMultiphysics.NeuralNetworkApplication as NeuralNetwork
from KratosMultiphysics.NeuralNetworkApplication.neural_network_process import NeuralNetworkProcess


from importlib import import_module
from tensorflow.keras import layers
import tensorflow.keras as keras

def Factory(settings):
    if not isinstance(settings, KM.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return AddLayerProcess(settings["parameters"])

class AddLayerProcess(NeuralNetworkProcess):

    def __init__(self, settings):
        super().__init__()
        """ The default constructor of the class

        Keyword arguments:
        self -- It signifies an instance of a class.
        model -- the container of the different model parts.
        settings -- Kratos parameters containing process settings.
        """

        default_settings = KM.Parameters("""{
            "help"                     : "This process adds a layer to a network.",
            "layer_type"               : "",
            "layer_parameters"         : {}
        }""")

        settings.ValidateAndAssignDefaults(default_settings)

        self.layer_type = settings["layer_type"].GetString()
        self.layer_parameters = settings["layer_parameters"]
        # Import layer classes
        module_name = "KratosMultiphysics.NeuralNetworkApplication." + self.layer_type
        class_module = import_module(module_name)
        self.layer_class = class_module.Factory(self.layer_parameters)

    def Add(self, model, hp = None):
        """ Processes to add layers to a network. """
        layer = self.layer_class.Build(hp = hp)
        if isinstance(model, keras.Sequential):
            model.add(layer)
            return model
        else:
            return layer(model)

