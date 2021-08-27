import KratosMultiphysics as KM
from KratosMultiphysics.NeuralNetworkApplication.neural_network_layer import NeuralNetworkLayerClass
import KratosMultiphysics.NeuralNetworkApplication.data_loading_utilities as DataLoadingUtilities

from tensorflow.keras import layers
import numpy as np
import h5py

def Factory(settings):
    if not isinstance(settings, KM.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return InputLayer(settings)

class InputLayer(NeuralNetworkLayerClass):
    """This class generates a base class for a Dense layer

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
            "data_input"               : "",
            "input_to_recurrent"       : false,
            "batch_size"               : "",
            "sparse"                   : false,
            "ragged"                   : false
        }""")

        settings.ValidateAndAssignDefaults(default_settings)

        # Get the input shape from a data input
        self.input_file = settings["data_input"].GetString()
        self.input_to_recurrent = settings["input_to_recurrent"].GetBool()

        # Input
        self.layer_name = settings["layer_name"].GetString()
        if not settings["batch_size"].GetString() == '':
            self.batch_size = settings["batch_size"].GetInt()
        else:
            self.batch_size = None
        self.sparse = settings["sparse"].GetBool()
        # self.tensor = settings["tensor"].GetString()
        self.ragged = settings["ragged"].GetBool()

        # When called by the add_layer_process, input the parameters in the keras function

    def Build(self, hp = None):
    
        self.input = DataLoadingUtilities.ImportDataFromFile(self.input_file, "InputData", lookback = self.input_to_recurrent).ExportAsArray()
        self.shape = self.input[0,:].shape
        self.layer = layers.Input(shape = self.shape, batch_size = self.batch_size, sparse = self.sparse,ragged=self.ragged,name=self.layer_name)
        
        return self.layer

