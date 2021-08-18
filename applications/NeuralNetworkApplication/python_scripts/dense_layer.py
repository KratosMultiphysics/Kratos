import KratosMultiphysics as KM
from KratosMultiphysics.NeuralNetworkApplication.neural_network_layer import NeuralNetworkLayerClass
from KratosMultiphysics.NeuralNetworkApplication.load_tunable_parameters_utility import LoadTunableParametersUtility
import KratosMultiphysics.NeuralNetworkApplication.data_loading_utilities as DataLoadingUtilities

from tensorflow.keras import layers
import tensorflow.keras as keras
import kerastuner as kt

def Factory(settings):
    if not isinstance(settings, KM.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return DenseLayer(settings)

class DenseLayer(NeuralNetworkLayerClass):
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
            "trainable"                : true,
            "output_layer"             : false,
            "data_output"              : "",
            "dtype"                    : "",
            "dynamic"                  : false,
            "units"                    : 1,
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
        self.output_layer = settings["output_layer"].GetBool()
        self.data_output = settings["data_output"].GetString()
        self.dtype = settings["dtype"].GetString()
        self.dynamic = settings["dynamic"].GetBool()
        if self.output_layer:
            self.data_output = settings["data_output"].GetString()
        else:
            self.units = settings["units"].GetInt()
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
        if self.output_layer:
            data = DataLoadingUtilities.ImportDataFromFile(self.data_output, "OutputData")
            self.units = data[0,:].size

        if self.tunable:
            setattr(self,self.tunable_variable,LoadTunableParametersUtility(self.tuning_parameters).Load(hp))
        
        self.layer = layers.Dense(self.units, activation=self.activation,use_bias=self.use_bias,
        kernel_initializer=self.kernel_initializer,bias_initializer=self.bias_initializer,kernel_regularizer=self.kernel_regularizer,
        bias_regularizer=self.bias_regularizer,activity_regularizer=self.activity_regularizer,
        trainable=self.trainable, dtype=self.dtype, name=self.layer_name, dynamic=self.dynamic)
        return self.layer



