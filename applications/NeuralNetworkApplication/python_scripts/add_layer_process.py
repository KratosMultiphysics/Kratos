import KratosMultiphysics as KM
import KratosMultiphysics.NeuralNetworkApplication as NeuralNetwork

from tensorflow.keras import layers

def Factory(settings, network):
    if not isinstance(settings, KM.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return AddLayerProcess(network, settings["Parameters"])

class AddLayerProcess(KM.Process):

    def __init__(self, network, settings):
        """ The default constructor of the class

        Keyword arguments:
        self -- It signifies an instance of a class.
        model -- the container of the different model parts.
        settings -- Kratos parameters containing process settings.
        """

        default_settings = KM.Parameters("""{
            "help"                     : "This process adds a layer to an initiatlized network.",
            "network_name"             : "",
            "layer_type"               : "",
            "layer_parameters"         : {}
        }""")

        settings.ValidateAndAssignDefaults(default_settings)

        self.network_name = settings["network_name"].GetString()
        self.layer_type = settings["layer_type"].GetString()
        self.layer_parameters = settings["layer_parameters"].GetParameters()

        # Load network

        # Switch with layer types maybe dynamic with the type_layer (e.g. dense_layer)

        # Save network
