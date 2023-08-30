import KratosMultiphysics as KM
import KratosMultiphysics.NeuralNetworkApplication.data_loading_utilities as DataLoadingUtilities
from KratosMultiphysics.NeuralNetworkApplication.neural_network_layer import NeuralNetworkLayerClass

from tensorflow.keras import layers
import tensorflow.keras as keras

def Factory(settings):
    if not isinstance(settings, KM.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return AddDimensionLayer(settings)

class AddDimensionLayer(NeuralNetworkLayerClass):
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
            "new_dimensions_sizes" : [1],
            "data_input" : "",
            "input_to_recurrent" : true
        }""")
        settings.ValidateAndAssignDefaults(default_settings)
        new_dimensions_sizes = settings['new_dimensions_sizes'].GetVector()
        self.new_size = tuple((-1,))
        for dimension in new_dimensions_sizes:
            self.new_size = self.new_size + tuple((int(dimension),))
        self.input_file = settings["data_input"].GetString()
        self.input_to_recurrent = settings["input_to_recurrent"].GetBool()


    def Build(self, hp = None):
        """ This method builds the layer when called by a process."""
        self.input = DataLoadingUtilities.ImportDataFromFile(self.input_file, "InputData", lookback = self.input_to_recurrent).ExportAsArray()
        self.shape = self.input[0,:].shape
        
        self.layer = layers.Reshape(self.shape[:-1]+self.new_size)
        return self.layer



