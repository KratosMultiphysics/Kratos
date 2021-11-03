import KratosMultiphysics as KM
from KratosMultiphysics.NeuralNetworkApplication.neural_network_layer import NeuralNetworkLayerClass
from KratosMultiphysics.NeuralNetworkApplication.load_tunable_parameters_utility import LoadTunableParametersUtility

from tensorflow.keras import layers
import tensorflow.keras as keras

def Factory(settings):
    if not isinstance(settings, KM.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return AveragePooling2DLayer(settings)

class AveragePooling2DLayer(NeuralNetworkLayerClass):
    """This class generates a base class for an MaxPooling1D layer.

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
            "pool_size_1"              : 2,
            "pool_size_2"              : 2,
            "strides_1"                : 0,
            "strides_2"                : 0,
            "padding"                  : "valid",
            "data_format"              : "channels_last",
            "tunable"                     : false,
            "tunable_variable"            : "",
            "tuning_parameters"           : {}
        }""")
        settings.ValidateAndAssignDefaults(default_settings)

        self.tunable = settings["tunable"].GetBool()
        self.layer_name = settings["layer_name"].GetString()

        if not settings["strides_1"].GetInt() == 0:
            strides_1 = settings["strides"].GetInt()
        elif not self.tunable:
            strides_1 = None
        
        if not settings["strides_2"].GetInt() == 0:
            strides_2 = settings["strides"].GetInt()
        elif not self.tunable:
            strides_2 = None

        if strides_1 == None or strides_2 == None:
            self.strides = None
        else:
            self.strides = tuple((strides_1, strides_2))
        self.pool_size = tuple((settings["pool_size_1"].GetInt(),settings["pool_size_2"].GetInt()))
        self.padding = settings["padding"].GetString()
        self.data_format = settings["data_format"].GetString()

        if self.tunable:
            self.tunable_variable = settings["tunable_variable"].GetString()
            self.tuning_parameters = settings["tuning_parameters"]

    def Build(self, hp = None):
        """ This method builds the layer when called by a process."""
        if self.tunable:
            setattr(self,self.tunable_variable,LoadTunableParametersUtility(self.tuning_parameters).Load(hp))
        self.layer = layers.AveragePooling2D(pool_size=self.pool_size, strides=self.strides, padding=self.padding, 
        data_format=self.data_format, name=self.layer_name)
        return self.layer



