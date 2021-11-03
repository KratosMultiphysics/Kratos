import KratosMultiphysics as KM
import KratosMultiphysics.NeuralNetworkApplication.data_loading_utilities as DataLoadingUtilities
from KratosMultiphysics.NeuralNetworkApplication.neural_network_layer import NeuralNetworkLayerClass

from tensorflow.keras import layers
import tensorflow.keras as keras

def Factory(settings):
    if not isinstance(settings, KM.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return ReshapeLayer(settings)

class ReshapeLayer(NeuralNetworkLayerClass):
    """This class generates a base class for a Reshape layer
    that adds a dimension to the input.

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
            "new_shape" : [1]
        }""")
        settings.ValidateAndAssignDefaults(default_settings)
        new_shape_sizes = settings['new_shape'].GetVector()
        self.new_shape = tuple()
        for dimension in new_shape_sizes:
            self.new_shape = self.new_shape + tuple((int(dimension),))



    def Build(self, hp = None):
        """ This method builds the layer when called by a process."""
        
        self.layer = layers.Reshape(self.new_shape)
        return self.layer



