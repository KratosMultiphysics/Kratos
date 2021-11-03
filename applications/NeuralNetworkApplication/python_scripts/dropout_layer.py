import numpy
import KratosMultiphysics as KM
from KratosMultiphysics.NeuralNetworkApplication.neural_network_layer import NeuralNetworkLayerClass

from tensorflow.keras import layers
import tensorflow.keras as keras

def Factory(settings):
    if not isinstance(settings, KM.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return DropoutLayer(settings)

class DropoutLayer(NeuralNetworkLayerClass):
    """This class generates a base class for a Batch Normalization layer

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
            "rate": 0.0,
            "noise_shape": [],
            "seed": "None"
        }""")
        settings.ValidateAndAssignDefaults(default_settings)

        self.rate = settings["rate"].GetDouble()
        if settings["noise_shape"].GetVector().Size() > 0:
            self.noise_shape = tuple()
            for dimension in settings["noise_shape"].GetVector():
                self.noise_shape = self.noise_shape + tuple((dimension,))
        else:
            self.noise_shape = None
        if settings["seed"].IsString():
            self.seed = None
        else:
            self.seed = settings["seed"].GetInt()


    def Build(self, hp = None):
        """ This method builds the layer when called by a process."""
        self.layer = layers.Dropout(self.rate, noise_shape=self.noise_shape, seed= self.seed)
        return self.layer



