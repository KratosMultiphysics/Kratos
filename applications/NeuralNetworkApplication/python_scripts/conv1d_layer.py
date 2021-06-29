import KratosMultiphysics as KM
from KratosMultiphysics.NeuralNetworkApplication.neural_network_layer import NeuralNetworkLayerClass
from KratosMultiphysics.NeuralNetworkApplication.load_tunable_parameters_utility import LoadTunableParametersUtility

from tensorflow.keras import layers
import tensorflow.keras as keras
import kerastuner as kt

def Factory(settings):
    if not isinstance(settings, KM.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return Conv1DLayer(settings)

class Conv1DLayer(NeuralNetworkLayerClass):
    """This class generates a base class for a Conv1D layer

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
            "trainable"                : true,
            "dtype"                    : "",
            "dynamic"                  : false,
            "filters"                  : 1,
            "kernel_size"              : 1,
            "strides"                  : 1,
            "padding"                  : "valid",
            "data_format"              : "channels_last",
            "dilation_rate"            : 1,
            "groups"                   : 1,
            "activation"               : "",
            "use_bias"                 : true,
            "kernel_initializer"       : "glorot_uniform",
            "bias_initializer"         : "zeros",
            "kernel_regularizer_l1"       : 0.0,
            "bias_regularizer_l1"         : 0.0,
            "activity_regularizer_l1"     : 0.0,
            "kernel_regularizer_l2"       : 0.0,
            "bias_regularizer_l2"         : 0.0,
            "activity_regularizer_l2"     : 0.0,
            "tunable"                     : false,
            "tunable_variable"            : "",
            "tuning_parameters"           : {}
        }""")
    # Constraints are currently not supported. Every supported constrain in keras should be implemented individually
        settings.ValidateAndAssignDefaults(default_settings)

        self.layer_name = settings["layer_name"].GetString()
        self.trainable = settings["trainable"].GetBool()
        self.dtype = settings["dtype"].GetString()
        self.dynamic = settings["dynamic"].GetBool()
        self.filters = settings["filters"].GetInt()
        self.kernel_size = settings["kernel_size"].GetInt()
        self.strides = settings["strides"].GetInt()
        self.padding = settings["padding"].GetString()
        self.data_format = settings["data_format"].GetString()
        self.dilation_rate = settings["dilation_rate"].GetInt()
        self.groups = settings["groups"].GetInt()
        if not settings["activation"].GetString() == "":
            self.activation = settings["activation"].GetString()
        else:
            self.activation = None
        self.use_bias = settings["use_bias"].GetBool()
        self.kernel_initializer = settings["kernel_initializer"].GetString()
        self.kernel_regularizer_l1 = settings["kernel_regularizer_l1"].GetDouble()
        self.kernel_regularizer_l2 = settings["kernel_regularizer_l2"].GetDouble()
        self.kernel_regularizer = keras.regularizers.L1L2(l1=self.kernel_regularizer_l1,l2=self.kernel_regularizer_l2)
        self.bias_initializer = settings["bias_initializer"].GetString()
        self.bias_regularizer_l1 = settings["bias_regularizer_l1"].GetDouble()
        self.bias_regularizer_l2 = settings["bias_regularizer_l2"].GetDouble()
        self.bias_regularizer = keras.regularizers.L1L2(l1=self.bias_regularizer_l1,l2=self.bias_regularizer_l2)
        self.activity_regularizer_l1 = settings["activity_regularizer_l1"].GetDouble()
        self.activity_regularizer_l2 = settings["activity_regularizer_l2"].GetDouble()
        self.activity_regularizer = keras.regularizers.L1L2(l1=self.activity_regularizer_l1,l2=self.activity_regularizer_l2)

        self.tunable = settings["tunable"].GetBool()
        
        if self.tunable:
            self.tunable_variable = settings["tunable_variable"].GetString()
            self.tuning_parameters = settings["tuning_parameters"]

    def Build(self, hp = None):
        """ This method builds the layer when called by a process."""
        if self.tunable:
            setattr(self,self.tunable_variable,LoadTunableParametersUtility(self.tuning_parameters).Load(hp))
        
        self.layer = layers.Conv1D(self.filters, self.kernel_size, strides=self.strides, padding=self.padding, 
        data_format=self.data_format, dilation_rate=self.dilation_rate, groups= self.groups, activation=self.activation,
        use_bias=self.use_bias, kernel_initializer=self.kernel_initializer,bias_initializer=self.bias_initializer,
        kernel_regularizer=self.kernel_regularizer, bias_regularizer=self.bias_regularizer,activity_regularizer=self.activity_regularizer,
        trainable=self.trainable, dtype=self.dtype, name=self.layer_name, dynamic=self.dynamic)
        return self.layer



