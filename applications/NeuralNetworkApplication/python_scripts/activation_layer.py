import KratosMultiphysics as KM
from KratosMultiphysics.NeuralNetworkApplication.neural_network_layer import NeuralNetworkLayerClass
from KratosMultiphysics.NeuralNetworkApplication.load_tunable_parameters_utility import LoadTunableParametersUtility

from tensorflow.keras import layers
import tensorflow.keras as keras

def Factory(settings):
    if not isinstance(settings, KM.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return ActivationLayer(settings)

class ActivationLayer(NeuralNetworkLayerClass):
    """This class generates a base class for an Activation layer.

    Public member variables:
    settings -- Kratos parameters containing process settings.
    """
    def __init__(self, settings):
        """ The default constructor of the class

        Keyword arguments:
        self -- It signifies an instance of a class.
        settings -- Kratos parameters containing layer settings.
        """

        default_settings = KM.Parameters("""{
            "layer_name"               : "",
            "activation"               : "",
            "tunable"                     : false,
            "tunable_variable"            : "",
            "tuning_parameters"           : {}
        }""")
        settings.ValidateAndAssignDefaults(default_settings)

        self.tunable = settings["tunable"].GetBool()
        self.layer_name = settings["layer_name"].GetString()
        if self.layer_name == "":
            self.layer_name = None
        if not settings["activation"].GetString() == "":
            self.activation = settings["activation"].GetString()
        elif not self.tunable:
            raise Exception("No activation function specified.")

        if self.tunable:
            self.tunable_variable = settings["tunable_variable"].GetString()
            self.tuning_parameters = settings["tuning_parameters"]

    def Build(self, hp = None):
        """ This method builds the layer when called by a process."""
        if self.tunable:
            setattr(self,self.tunable_variable,LoadTunableParametersUtility(self.tuning_parameters).Load(hp))
        self.layer = layers.Activation(self.activation, name=self.layer_name)
        return self.layer



