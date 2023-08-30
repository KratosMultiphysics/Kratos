import KratosMultiphysics as KM
from KratosMultiphysics.NeuralNetworkApplication.neural_network_layer import NeuralNetworkLayerClass

from tensorflow.keras import layers
import tensorflow.keras as keras

def Factory(settings):
    if not isinstance(settings, KM.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return RepeatVectorLayer(settings)

class RepeatVectorLayer(NeuralNetworkLayerClass):
    """This class generates a base class for a RepeatVector layer

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
            "n" : 1
        }""")
        settings.ValidateAndAssignDefaults(default_settings)
        self.n = settings["n"].GetInt()


    def Build(self, hp = None):
        """ This method builds the layer when called by a process."""
        self.layer = layers.RepeatVector(self.n)
        return self.layer



