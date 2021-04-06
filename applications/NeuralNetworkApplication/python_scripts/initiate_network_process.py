import KratosMultiphysics as KM
import KratosMultiphysics.NeuralNetworkApplication as NeuralNetwork

from tensorflow.keras import layers

def Factory(settings, network):
    if not isinstance(settings, KM.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return InitiateNetworkProcess(network, settings["Parameters"])

class InitiateNetworkProcess(KM.Process):

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
            "network_type"             : "Sequential"
            "initial_layer_type"       : "",
            "initial_layer_parameters" : {}
        }""")

        settings.ValidateAndAssignDefaults(default_settings)

        self.network_name = settings["network_name"].GetString()
        self.network_type = settings["network_type"].GetString()
        self.initial_layer_type = settings["initial_layer_type"].GetString()
        self.initial_layer_parameters = settings["initial_layer_parameters"].GetParameters()

        # Create network in path

        # Initiate the network with the given layer and type (Sequential or functional)

        # Save network